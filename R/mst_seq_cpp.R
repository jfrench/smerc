# library(Rcpp)
#
# stat_pois <- '
# Rcpp::NumericVector stat_poisson_cp(Rcpp::NumericVector yin,
#                                double ty,
#                                Rcpp::NumericVector ein,
#                                Rcpp::NumericVector eout,
#                                unsigned int min_cases,
#                                Rcpp::NumericVector popin,
#                                double max_pop) {
#   unsigned int yin_length = yin.length();
#   double lrin, lrout;
#   Rcpp::NumericVector tall(yin_length, 0);
#   Rcpp::NumericVector yout(yin_length, 0);
#
#   // determine if there will be any problematic statistics
#   for (unsigned int i = 0; i < yin_length; i++) {
#     // compute statistic for good locations
#     // yin > 0 and yin/ein > yout/ein
#     if (yin[i] >= min_cases & popin[i] <= max_pop) {
#       yout[i] = ty - yin[i];
#       lrin = log(yin[i]) - log(ein[i]);
#       lrout = log(yout[i]) - log(eout[i]);
#       if (lrin > lrout) {
#         tall[i] = yin[i] * lrin + yout[i] * lrout;
#       }
#     }
#   }
#
#   return tall;
# }
# '
#
# cppFunction(stat_pois)
#
# smerc::stat.poisson(
#   yin = 106, yout = 552 - 106,
#   ein = 62.13, eout = 552 - 62.13
# )
#
# stat_poisson_cp(106,
#                 552,
#                 62.13,
#                 552 - 62.13,
#                 2,
#                 15000,
#                 20000)
#
#
# mst_seq <- 'Rcpp::List(int region,
#                       Rcpp::IntegerVector neighbors,
#                       Rcpp::NumericVector cases,
#                       Rcpp::NumericVector pop,
#                       Rcpp::NumericVector ex,
#                       Rcpp::IntegerMatrix w,
#                       double ty,
#                       double max_pop,
#                       unsigned int type,
#                       unsigned int nlinks,
#                       bool early
#                       unsigned int min_cases) {
#   // initialize objects
#   unsigned int nneighbs = neighbors.size();
#   loglikrat = Rcpp::NumericVector(nneighbs);
#   yin = Rcpp::NumericVector(nneighbs);
#   ein = Rcpp::NumericVector(nneighbs);
#   popin = Rcpp::NumericVector(nneighbs);
#
#   // first step
#   yin[0] = cases[region];
#   ein[0] = ex[region];
#   popin[0] = pop[region];
#   loglikrat = stat_poisson_cp(yin, ty, ein, eout, min_cases, popin, max_pop);
#
#   // loglikrat[1] <- scan.stat(yin[1], ein[1], ty - ein[1], ty)
#   // uz <- max_neighbors <- vector("list", length(neighbors))
#   // uz[[1]] <- region
#   // max_neighbors[[1]] <- setdiff(neighbors, region)
#
#   return List::create(
#     Named("loglikrat") = loglikrat,
#     Named("cases") = yin,
#     Named("expected") = ein,
#     Named("population") = popin
#   );
# }
# '
#
# cppFunction(code = mst_seq, includes = stat_pois)
#
#
#
# start, neighbors, cases, pop, w, ex, ty,
# max_pop, type = "maxonly", nlinks = "one",
# early = FALSE) {
#   loglikrat <- yin <- ein <- popin <- numeric(length(neighbors))
#   region <- start
#   yin[1] <- cases[region]
#   ein[1] <- ex[region]
#   popin[1] <- pop[region]
#   loglikrat[1] <- scan.stat(yin[1], ein[1], ty - ein[1], ty)
#   uz <- max_neighbors <- vector("list", length(neighbors))
#   uz[[1]] <- region
#   max_neighbors[[1]] <- setdiff(neighbors, region)
#
# # body5 <-
# # '
# # arma::mat cov_spBayes(arma::mat D, int sp_type, double sigmasq,
# # 	double phi, double nu, double ev, double fv)
# # {
# # 	int nr = D.n_rows;
# # 	int nc = D.n_cols;
# #
# # 	arma::mat V = arma::zeros(nc, nc);
# #
# # 	for(int i = 0; i < nr; i++)
# # 	{
# # 		for(int j = 0; j < nc; j++)
# # 		{
# # 			if(D(i, j) == 0)
# # 			{
# # 				V(i, j) = sigmasq + ev + fv;
# # 			}
# # 			else
# # 			{
# # 				if(sp_type == 1)
# # 				{
# # 					V(i, j) = sigmasq * exp(-D(i, j) * phi);
# # 				}
# # 				else if(sp_type == 2)
# # 				{
# # 					V(i, j) = sigmasq * exp(-pow(D(i, j) * phi, 2));
# # 				}
# # 				else if(sp_type == 3)
# # 				{
# # 					V(i, j) = sigmasq * pow(D(i,j)*phi, nu)/(pow(2, nu-1)*R::gammafn(nu))*R::bessel_k(D(i,j)*phi, nu, 1.0);
# # 				}
# # 				else
# # 				{
# # 					if(D(i, j) <= 1.0/phi)
# # 					{
# # 						V(i, j) = sigmasq * (1.0 - 1.5*phi*D(i,j) + 0.5*pow(phi*D(i,j),3));
# # 					}
# #
# # 				}
# # 			}
# # 		}
# # 	}
# # 	return V;
# # }
# # '
# #
# # cppFunction(body5, depends = "RcppArmadillo")
