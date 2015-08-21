#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#define TOL 0.00000001

using namespace Rcpp;

// Currently Rcpp have no conversion to cube, so we need this function
arma::cube array_to_cube(NumericVector v, int x, int y, int z) {
    arma::cube cube_array(v.begin(), x, y, z, false);
    return cube_array;
}

arma::cube get_E(int n_states, arma::mat qvec, arma::mat qvec_inv) {
    arma::cube E(n_states, n_states, n_states);
    arma::mat E_i(n_states, n_states, arma::fill::zeros);
    int i;

    for (i=0; i<n_states; i++) {
        E_i(i,i) = 1;
        E.slice(i) = qvec * E_i * qvec_inv;
        E_i(i,i) = 0;
    }

    return E;
}

arma::cube get_Si(int n_states, arma::cube E, arma::mat m) {
    arma::cube Si(n_states, n_states, n_states);
    for (int i=0; i<n_states; i++) {
        Si.slice(i) = E.slice(i) * m;
    }
    return Si;
}

double build_Iij(double t, arma::vec d, int i, int j) {
    double diff = d[i] - d[j];
    if (abs(diff) < TOL) {
        return t * (exp(d[i]*t));
    } else {
        return (exp(d[i]*t) - exp(d[j]*t)) / (diff);
    }
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::cube transition_probabilities(List Q_eigen, arma::vec edges, int omp) {

    arma::vec q_val = Q_eigen["values"];
    arma::mat q_vec_inv = Q_eigen["vectors_inv"];
    arma::mat q_vec = Q_eigen["vectors"];

    int n = q_val.size();
    arma::cube P(n, n, edges.size());

    omp_set_num_threads(omp);
    #pragma omp parallel default(shared)
    {
        arma::mat aux(n, n);
        int i, j;

        #pragma omp for nowait
        for (i=0; i<edges.size(); i++) {
            for (j=0; j<n; j++) {
                aux.col(j) = q_vec.col(j) * exp(q_val(j) * edges(i));
            }
            P.slice(i) = aux * q_vec_inv;
        }
    }

    return (P);
}

// The expected number of labelled markov transitions and expected
// markov rewards
// TODO: find a better name for this function =P
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List func_H(arma::mat Q, List Q_eigen, List tree, List tree_extra, int omp) {
    int n_edges = tree_extra["n_edges"];
    int n_states = tree_extra["n_states"];
    arma::cube lmt(n_states, n_states, n_edges, arma::fill::zeros);
    arma::cube emr(n_states, n_states, n_edges, arma::fill::zeros);
    arma::vec d = Q_eigen["values"];
    arma::mat vectors = Q_eigen["vectors"];
    arma::mat vectors_inv = Q_eigen["vectors_inv"];
    arma::mat emr_diag = arma::diagmat(as<arma::vec>(tree_extra["rewards"]));
    arma::vec edges = tree["edge.length"];
    arma::mat QL = Q;
    QL.diag().zeros();
    List ret;

    arma::cube E = get_E(n_states, vectors, vectors_inv);
    arma::cube s_lmt = get_Si(n_states, E, QL);
    arma::cube s_emr = get_Si(n_states, E, emr_diag);

    omp_set_num_threads(omp);
    #pragma omp parallel default(shared)
    {
        int b, i, j;
        double Iij, edge;

        #pragma omp for nowait
        for (b=0; b<n_edges; b++) {
            edge = edges[b];
            for (i=0; i<n_states; i++) {
                for (j=0; j<n_states; j++) {
                    Iij = build_Iij(edge, d, i, j);
                    lmt.slice(b) += s_lmt.slice(i) * E.slice(j) * Iij;
                    emr.slice(b) += s_emr.slice(i) * E.slice(j) * Iij;
                }
            }
        }
    }

    ret["lmt"] = lmt;
    ret["emr"] = emr;
    return ret;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List posterior_restricted_moment(List tree, List tree_extra, List map, int omp) {
    int n_edges = tree_extra["n_edges"];
    int n_states = tree_extra["n_states"];
    arma::mat prm_lmt(n_edges, n_states, arma::fill::zeros);
    arma::mat prm_emr(n_edges, n_states, arma::fill::zeros);

    List h = map["h"];
    List fl = map["fl"];

    arma::cube emr = array_to_cube(h["emr"], n_states, n_states, n_edges);
    arma::cube lmt = array_to_cube(h["lmt"], n_states, n_states, n_edges);
    arma::mat f = fl["F"];
    arma::mat g = fl["G"];
    arma::mat s = fl["S"];
    arma::mat edges = tree["edge"];

    List ret;

    omp_set_num_threads(omp);
    #pragma omp parallel default(shared)
    {
        int i, j, p, c, b;
        arma::mat lmt_i(n_states, n_states);
        arma::mat emr_i(n_states, n_states);
        arma::rowvec gs(n_states);

        #pragma omp for nowait
        for (i=0; i<n_edges; i++) {
            p = edges(i,0)-1; // nodes start at 1 in R..
            c = edges(i,1)-1; // nodes start at 1 in R..
            b = (i%2==0? edges(i+1,1) : edges(i-1,1)) - 1;

            lmt_i = lmt.slice(i);
            emr_i = emr.slice(i);

            for (j=0; j<n_states; j++) {
                gs = g(p,j) * s(b,j) * f.row(c);

                prm_lmt(i,j) += (sum(gs % lmt_i.row(j)));
                prm_emr(i,j) += (sum(gs % emr_i.row(j)));
            }
        }
    }

    ret["lmt"] = prm_lmt;
    ret["emr"] = prm_emr;
    return ret;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List fractional_likelihoods(List tree, List tree_extra, arma::mat q
                            , List q_eigen, arma::rowvec prior
                            , NumericVector trans_prob, double tol)
{

    int n_edges = tree_extra["n_edges"];
    int n_states = tree_extra["n_states"];
    int n_nodes = tree_extra["n_nodes"];
    int n_tips = tree_extra["n_tips"];
    arma::cube tp = array_to_cube(trans_prob, n_states, n_states, n_edges);

    arma::mat f(n_nodes, n_states, arma::fill::zeros);
    arma::mat g(n_nodes, n_states, arma::fill::zeros);
    arma::mat s(n_nodes, n_states, arma::fill::zeros);
    arma::mat edge = tree["edge"];
    arma::mat states = tree_extra["states"];
    arma::mat t_right(n_states, n_states);
    arma::mat t_left(n_states, n_states);
    double likelihood;
    int i, j, e, p, right, left, root_node;
    List ret;

    // Initialize f
    for (i=0; i<n_tips; i++) {
        f.row(i) = states.row(i);
    }

    for (e=0; e<n_edges; e+=2) {
        p = edge(e,0) - 1; // nodes start at 1 in R...
        right = edge(e,1) - 1;
        left = edge(e+1,1) - 1;

        t_right = tp.slice(e);
        t_left = tp.slice(e+1);

        for (i=0; i<n_states; i++) {
            for (j=0; j<n_states; j++) {
                s(right,i) += f(right,j) * t_right(i,j);
                s(left,i)  += f(left,j) * t_left(i,j);
            }
            f(p,i) = s(right,i) * s(left,i);
        }
    }

    root_node = edge(n_edges-1,0)-1;
    likelihood = sum(f.row(root_node) % prior);
    g.row(root_node) = prior;

    for (e=n_edges-1; e>=1; e-=2) {
        p = edge(e,0) - 1;
        left = edge(e,1) - 1;
        right = edge(e-1,1) - 1;

        t_left = tp.slice(e);
        t_right = tp.slice(e-1);

        for (i=0; i<n_states; i++) {
            for (j=0; j<n_states; j++) {
                g(left,i)  += g(p,j) * s(right,j) * t_left(i,j);
                g(right,i) += g(p,j) * s(left,j) * t_right(i,j);
            }
        }
    }

    ret["F"] = f;
    ret["G"] = g;
    ret["S"] = s;
    ret["L"] = likelihood;

    return(ret);
}
