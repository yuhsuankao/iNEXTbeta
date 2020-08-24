#' \code{iNEXTbeta}: compute the rarefaction and extrapolation of beta diveristy.
#' @param x a list consist of N data.frame/matrix describing species-by-assemblage abundance, if the data_type is "abundance".
#'          a list consist of N lists and each element of the N list is a list consist of T data.frame/matrix describing raw incidence of every assemblage, if the data_type is "incidence_raw".
#'          Note that the species in each element must exactly match including species order.
#' @param coverage_expected a numeric vector of sample coverage for which diversity estimates will be computed. Each number should between 0 and 1.
#' @param data_type data type of input data: individual-based abundance data (\code{data_type = "abundance"}) or sampling-unit-based incidence data (\code{data_type = "incidence_raw"})
#' @param level the level of diversity to be computed.
#' @param nboot an integer specifying the number of replications for se estimating. Setting 0 to ignore se estimating.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param phy_tree a phylo object describing the phylogenetic tree in Newick format for all observed species in the pooled assemblage. used if \code{level="phylogenetic"}.
#' @param reftime a positive value specifying the reference times for diversity computation. If \code{NULL},
#'                then \code{reftime} is set to be the tree depth of the phylogenetic tree, which is spanned by all the observed species in
#'                the pooled assemblage. Default is \code{NULL}. used if \code{level="phylogenetic"}.
#' @param distance_matrix a matrix describing the functional distance for all observed species in the pooled assemblage. used if \code{level="functional"}.
#' @param tau_type the functional diversity should be computed at given tau value \code{tau_type="single"} or for AUC \code{tau_type="AUC"}.
#' @param tau a positvie number specifying the level of threshold distinctiveness. used if \code{level="functional" and tau_type="single"}.
#' @param cut_number a integer specifying the number of cut when computing the AUC of functional diversity. used if \code{level="functional" and tau_type="AUC"}.

#' @return a list consists of N elements, where N is the number of region. Each element is a list containing 7 tables of gamma, alpha and beta diversities,
#' and 4 types of dissimialrities (1-CqN, 1-UqN, 1-SqN, 1-VqN), respectively.
#' @import tidyverse
#' @import magrittr
#' @import Rcpp
#' @import future.apply
#' @import readxl
#' @import abind
#' @import iNEXT
#' @import FunD
#' @import PhD
#' @import ade4
#' @import phytools
#' @import phyclust
#' @import chaoUtility
#' @import tidytree
#' @export

iNEXTbeta = function(x, coverage_expected, data_type = c('abundance', 'incidence_raw'), level = c('taxonomic', 'phylogenetic', 'functional'),
                     nboot = 20, conf = 0.95,
                     phy_tree = NULL, reftime = NULL,
                     distance_matrix = NULL, tau_type = c('single', 'AUC'), tau = NULL, cut_number = NULL){

  by = 'coverage'

  if(data_type=='abundance'){

    if( class(x)=="data.frame" | class(x)=="matrix" ) x = list(Region_1 = x)

    if(class(x)== "list"){
      if(is.null(names(x))) region_names = paste0("Region_", 1:length(x)) else region_names = names(x)
      Ns = sapply(x, ncol)
      data_list = x
    }

  }

  if(data_type=='incidence_raw'){

    if(is.null(names(x))) region_names = paste0("Region_", 1:length(x)) else region_names = names(x)
    Ns = sapply(x, length)
    data_list = x

  }


  if(is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)

  if(level=='phylogenetic'){

    if (data_type=='abundance') {

      pool.data = do.call(cbind, data_list) %>% rowSums
      pool.data = do.call(cbind, data_list) %>% rowSums
    }

    if (data_type=='incidence_raw') {

      pool.data = do.call(cbind, data_list[[1]]) %>% rowSums
      pool.data = do.call(cbind, data_list[[1]]) %>% rowSums

    }

    pool.name = names(pool.data[pool.data>0])
    tip = phy_tree$tip.label[-match(pool.name, phy_tree$tip.label)]
    mytree = drop.tip(phy_tree, tip)
    H_max = get.rooted.tree.height(mytree)

    if(is.null(reftime)) { reft = H_max
    } else if (reftime <= 0) { stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
    } else { reft = reftime }

  }

  for_each_region = function(data, region_name, N){

    #data
    if (data_type=='abundance') {

      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector

      ref_gamma = iNEXT:::Chat.Ind(data_gamma, n)
      if (by=='size') ref_alpha = ref_gamma
      if (by=='coverage') ref_alpha = iNEXT:::Chat.Ind(data_alpha, n)
      # ref_alpha_max = iNEXT:::Chat.Ind(data_alpha, n*2)

      coverage_expected = coverage_expected[coverage_expected<ref_gamma]
      coverage_expected = c(coverage_expected, ref_gamma, ref_alpha) %>% sort %>% unique

      m_gamma = sapply(coverage_expected, function(i) coverage_to_size(data_gamma, i, data_type='abundance'))
      if (by=='size') m_alpha = m_gamma
      if (by=='coverage') m_alpha = sapply(coverage_expected, function(i) coverage_to_size(data_alpha, i, data_type='abundance'))

    }

    if (data_type=='incidence_raw') {

      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all regions must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)

      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))

      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)

      # data_gamma_freq = data_gamma_freq[data_gamma_freq>0]
      # data_alpha_freq = data_alpha_freq[data_alpha_freq>0]

      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame

      ref_gamma = iNEXT:::Chat.Sam(data_gamma_freq, n)
      if (by=='size') ref_alpha = ref_gamma
      if (by=='coverage') ref_alpha = iNEXT:::Chat.Sam(data_alpha_freq, n)
      # ref_alpha_max = iNEXT:::Chat.Sam(data_alpha_freq, n*2)

      coverage_expected = coverage_expected[coverage_expected<ref_gamma]
      coverage_expected = c(coverage_expected, ref_gamma, ref_alpha) %>% sort %>% unique

      m_gamma = sapply(coverage_expected, function(i) coverage_to_size(data_gamma_freq, i, data_type='incidence_freq'))
      if (by=='size') m_alpha = m_gamma
      if (by=='coverage') m_alpha = sapply(coverage_expected, function(i) coverage_to_size(data_alpha_freq, i, data_type='incidence_raw'))

    }



    if (level=='taxonomic') {

      if (data_type=='abundance') {

        gamma = lapply(1:length(coverage_expected), function(i){
          iNEXT::estimateD(as.numeric(data_gamma), datatype = "abundance", base = "coverage", level = coverage_expected[i], conf = NULL)
        }) %>% do.call(rbind,.)

        if (by=='size') {

          alpha = lapply(1:length(coverage_expected), function(i){
            iNEXT::estimateD(as.numeric(data_alpha), datatype = "abundance", base = "size", level = m_alpha[i], conf = NULL)
          }) %>% do.call(rbind,.)

        }
        if (by=='coverage') {

          alpha = lapply(1:length(coverage_expected), function(i){
            iNEXT::estimateD(as.numeric(data_alpha), datatype = "abundance", base = "coverage", level = coverage_expected[i], conf = NULL)
          }) %>% do.call(rbind,.)

        }

      }

      if (data_type=='incidence_raw') {

        gamma = lapply(1:length(coverage_expected), function(i){
          iNEXT::estimateD(as.numeric(data_gamma_freq), datatype = "incidence_freq", base = "coverage", level = coverage_expected[i], conf = NULL)
        }) %>% do.call(rbind,.)

        if (by=='size') {

          alpha = lapply(1:length(coverage_expected), function(i){
            iNEXT::estimateD(data_alpha_freq, datatype = "incidence_freq", base = "size", level = m_alpha[i], conf = NULL)
          }) %>% do.call(rbind,.)

        }
        if (by=='coverage') {

          alpha = lapply(1:length(coverage_expected), function(i){
            iNEXT::estimateD(data_alpha_freq, datatype = "incidence_freq", base = "coverage", level = coverage_expected[i], conf = NULL)
          }) %>% do.call(rbind,.)

        }



      }

      gamma = (cbind(coverage_expected, gamma[,-2]) %>%
                 mutate(Method = ifelse(coverage_expected<ref_gamma, 'Interpolated', 'Ref_gamma')) %>%
                 gather(Order, Estimate, 4:6))[,c(6,5,4,1,3,2)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      for (i in 0:2) gamma$Order[gamma$Order==paste0('q = ', i)] = i
      gamma$Order = as.numeric(gamma$Order)

      # if (max_alpha_coverage==T) under_max_alpha = !((gamma$Order==0) & (gamma$Coverage_expected>ref_alpha_max)) else under_max_alpha = gamma$Coverage_expected>0
      # gamma = gamma[under_max_alpha,]



      alpha = (cbind(coverage_expected, alpha[,-2]) %>%
                 mutate(Method = ifelse(coverage_expected<ref_alpha, 'Interpolated',
                                        ifelse(coverage_expected==ref_alpha, 'Ref_alpha',
                                               ifelse(coverage_expected==ref_gamma, 'Ref_gamma', 'Extrapolated')))) %>%
                 gather(Order, Estimate, 4:6))[,c(6,5,4,1,3,2)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      alpha$Estimate = alpha$Estimate / N

      for (i in 0:2) alpha$Order[alpha$Order==paste0('q = ', i)] = i
      alpha$Order = as.numeric(alpha$Order)

      # alpha = alpha[under_max_alpha,]



      beta = alpha
      beta$Estimate = gamma$Estimate/alpha$Estimate

      C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
      U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
      V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
      S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))

      if(nboot>1){

        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"coverage_expected","N",'under_max_alpha',
        #                     'data_type', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))

        # plan(sequential)
        # plan(multiprocess)

        # se = parSapply(cl, 1:nboot, function(i){

        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){

          if (data_type=='abundance') {

            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))

            bootstrap_data_gamma = rowSums(bootstrap_sample)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
            bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]

            m_gamma = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_gamma, i, data_type='abundance'))

            gamma = lapply(1:length(coverage_expected), function(i){
              iNEXT::estimateD(as.numeric(bootstrap_data_gamma), datatype = "abundance", base = "coverage", level = coverage_expected[i], conf = NULL)
            }) %>% do.call(rbind,.)

            if (by=='size') {

              alpha = lapply(1:length(coverage_expected), function(i){
                iNEXT::estimateD(as.numeric(bootstrap_data_alpha), datatype = "abundance", base = "size", level = m_gamma[i], conf = NULL)
              }) %>% do.call(rbind,.)

            }
            if (by=='coverage') {

              alpha = lapply(1:length(coverage_expected), function(i){
                iNEXT::estimateD(as.numeric(bootstrap_data_alpha), datatype = "abundance", base = "coverage", level = coverage_expected[i], conf = NULL)
              }) %>% do.call(rbind,.)

            }



          }

          if (data_type=='incidence_raw') {

            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')

            raw = lapply(1:ncol(bootstrap_population), function(j){

              lapply(1:nrow(bootstrap_population), function(i) rbinom(n=n, size=1, prob=bootstrap_population[i,j])) %>% do.call(rbind,.)

            })

            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))

            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)

            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq>0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq>0]

            m_gamma = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_gamma_freq, i, data_type='incidence_freq'))

            gamma = lapply(1:length(coverage_expected), function(i){
              iNEXT::estimateD(bootstrap_data_gamma_freq, datatype = "incidence_freq", base = "coverage", level = coverage_expected[i], conf = NULL)
            }) %>% do.call(rbind,.)

            if (by=='size') {

              alpha = lapply(1:length(coverage_expected), function(i){
                iNEXT::estimateD(bootstrap_data_alpha_freq, datatype = "incidence_freq", base = "size", level = m_gamma[i], conf = NULL)
              }) %>% do.call(rbind,.)

            }
            if (by=='coverage') {

              alpha = lapply(1:length(coverage_expected), function(i){
                iNEXT::estimateD(bootstrap_data_alpha_freq, datatype = "incidence_freq", base = "coverage", level = coverage_expected[i], conf = NULL)
              }) %>% do.call(rbind,.)

            }



          }



          gamma = (gamma[,4:6] %>% gather(Order, Estimate))$Estimate#[under_max_alpha]

          alpha = (alpha[,4:6] %>% gather(Order, Estimate))$Estimate#[under_max_alpha]
          alpha = alpha / N

          beta = gamma/alpha

          order = rep(c(0,1,2), each=length(coverage_expected))#[under_max_alpha]

          beta = data.frame(Estimate=beta, order)

          C = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(1-order) - 1)/(N^(1-order)-1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(order-1) - 1)/(N^(order-1)-1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate-1)/(N-1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1)))$Estimate

          beta = beta$Estimate

          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along=3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        se = matrix(0, ncol = 7, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)

      }

    }

    if (level=='phylogenetic') {

      if (data_type=='abundance') {

        aL = phyBranchAL_Abu(phylo = phy_tree, data = data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

        gamma = PhD:::PhD.m.est(aL=aL_table_gamma, m=m_gamma, Q=c(0,1,2), datatype='abundance', nt=n) %>% t %>% as.data.frame %>%
          set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
          mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Ind(data_gamma, m_gamma), 3), Size=rep(m_gamma, 3))


        aL_table_alpha = c()

        for (i in 1:N){

          x = data[data[,i]>0,i]
          names(x) = rownames(data)[data[,i]>0]

          aL = phyBranchAL_Abu(phylo = phy_tree, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }


        qPDm = PhD:::PhD.m.est(aL=aL_table_alpha, m=m_alpha, Q=c(0,1,2), datatype='abundance', nt=n)
        qPDm = qPDm/N
        alpha = qPDm %>% t %>% as.data.frame %>%
          set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
          mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Ind(data_alpha, m_alpha), 3), Size=rep(m_alpha, 3))

      }

      if (data_type=='incidence_raw') {

        aL = phyBranchAL_Inc(phylo=phy_tree, data=as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

        gamma = PhD:::PhD.m.est(aL = aL_table_gamma, m = m_gamma, Q = c(0,1,2), datatype = 'incidence_raw', nt = n) %>% t %>% as.data.frame %>%
          set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
          mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Sam(data_gamma_freq, m_gamma), 3), Size=rep(m_gamma, 3))

        aL_table_alpha = c()

        for (i in 1:N){

          x = data[[i]]

          aL = phyBranchAL_Inc(phylo = phy_tree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }

        alpha = (PhD:::PhD.m.est(aL = aL_table_alpha, m = m_alpha, Q = c(0,1,2), datatype = 'incidence_raw', nt = n)/N) %>% t %>% as.data.frame %>%
          set_colnames(c(0,1,2)) %>% gather(Order, Estimate) %>%
          mutate(Coverage_expected=rep(coverage_expected, 3), Coverage_real=rep(iNEXT:::Chat.Sam(data_alpha_freq, m_alpha), 3), Size=rep(m_alpha, 3))


      }

      gamma = (gamma %>%
                 mutate(Method = ifelse(coverage_expected<ref_gamma, 'Interpolated', 'Ref_gamma')))[,c(2,1,6,3,4,5)] %>%
        set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      # if (max_alpha_coverage==T) under_max_alpha = !((gamma$Order==0) & (gamma$Coverage_expected>ref_alpha_max)) else under_max_alpha = gamma$Coverage_expected>0
      # gamma = gamma[under_max_alpha,]
      gamma$Order = as.numeric(gamma$Order)



      alpha = (alpha %>%
                 mutate(Method = ifelse(coverage_expected<ref_alpha, 'Interpolated',
                                        ifelse(coverage_expected==ref_alpha, 'Ref_alpha',
                                               ifelse(coverage_expected==ref_gamma, 'Ref_gamma', 'Extrapolated')))))[,c(2,1,6,3,4,5)] %>%
        set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      # alpha = alpha[under_max_alpha,]
      alpha$Order = as.numeric(alpha$Order)

      beta = alpha
      beta$Estimate = gamma$Estimate/alpha$Estimate

      C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
      U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
      V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
      S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))

      if(nboot>1){

        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"coverage_expected","N",'under_max_alpha',
        #                     'data_type', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))

        # plan(sequential)
        # plan(multiprocess)

        # se = parSapply(cl, 1:nboot, function(i){

        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){

          if (data_type=='abundance') {

            tree_bt = phy_tree

            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol=ncol(data))

            if ( nrow(p_bt) > nrow(data) & sum(unseen_p)>0 ){

              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data), unseen_name)

              bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              x_bt = bootstrap_sample

              rownames(x_bt) = rownames(p_bt)

              if ( sum(x_bt[-(1:nrow(data)),])>0 ){

                g0_hat = apply(data, 2, function(x){

                  n = sum(x)
                  f1 = sum(x==1)
                  f2 = sum(x==2)

                  aL = phyBranchAL_Abu(phylo = phy_tree, data = x, rootExtend = T, refT = reft)

                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  g1 = aL$branch.length[aL$branch.abun==1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun==2] %>% sum
                  g0_hat = ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
                  g0_hat

                })

                te = (x_bt[1:nrow(data),]*(data==0))>0
                used_length = sapply(1:ncol(data), function(i) {

                  if (sum(te[,i])==0) return(0) else {

                    phyBranchAL_Abu(phylo = phy_tree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i]==TRUE))) %>% select(branch.length) %>% sum

                  }

                })

                g0_hat = g0_hat-used_length
                g0_hat[g0_hat<0] = 0

                unseen_sample = x_bt[-(1:nrow(data)),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol=ncol(x_bt), byrow=T)

                L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i]>0)>0) (g0_hat[i] / nrow(unseen)) else 0 )

                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow=T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample)==0)] = 0

                for (i in 1:length(L0_hat)){

                  tip = list(edge=matrix(c(2,1),1,2),
                             tip.label=unseen_name[i],
                             edge.length=L0_hat[i],
                             Nnode=1)
                  class(tip) = "phylo"

                  tree_bt = tree_bt + tip

                }

              } else {

                x_bt = x_bt[1:nrow(data),]
                p_bt = p_bt[1:nrow(data),]

              }

            } else {

              p_bt = p_bt[1:nrow(data),]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              rownames(x_bt) = rownames(data)

            }

            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]

            m_gamma = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_gamma, i, data_type='abundance'))
            m_alpha = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_alpha, i, data_type='abundance'))

            aL = phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

            gamma = as.vector(PhD:::PhD.m.est(aL=aL_table_gamma, m=m_gamma, Q=c(0,1,2), datatype='abundance', nt=n) %>% t)


            aL_table_alpha = c()

            for (i in 1:N){

              # x = x_bt[x_bt[,i]>0,i]
              # names(x) = rownames(p_bt)[x_bt[,i]>0]

              x = x_bt[,i]
              names(x) = rownames(p_bt)
              x = x[x_bt[,i]>0]

              aL = phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

              aL_table_alpha = rbind(aL_table_alpha, aL_table)

            }

            alpha = as.vector((PhD:::PhD.m.est(aL=aL_table_alpha, m=m_alpha, Q=c(0,1,2), datatype='abundance', nt=n)/N) %>% t)

          }

          if (data_type=='incidence_raw') {

            tree_bt = phy_tree

            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])),] %>% matrix(ncol=N)

            if ( nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p)>0 ){

              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data[[1]])),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data[[1]]), unseen_name)

              raw = lapply(1:ncol(p_bt), function(j){

                lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)

              })

              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)

              if ( lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])),]) %>% do.call(sum,.)>0 ){

                R0_hat = sapply(data, function(x){

                  nT = ncol(x)
                  Q1 = sum(rowSums(x)==1)
                  Q2 = sum(rowSums(x)==2)

                  aL = phyBranchAL_Inc(phylo = phy_tree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)

                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  R1 = aL$branch.length[aL$branch.abun==1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun==2] %>% sum
                  R0_hat = ifelse( R2>((R1*Q2)/(2*Q1)) , ((nT-1)/nT)*(R1^2/(2*R2)) , ((nT-1)/nT)*(R1*(Q1-1)/(2*(Q2+1))) )
                  R0_hat

                })

                te = (sapply(raw, rowSums)[1:nrow(data[[1]]),]*(sapply(data, rowSums)==0))>0
                used_length = sapply(1:N, function(i) {

                  if (sum(te[,i])==0) return(0) else {

                    phyBranchAL_Inc(phylo = phy_tree, data = raw[[i]][1:nrow(data[[1]]),], datatype = "incidence_raw", rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i]==TRUE))) %>% select(branch.length) %>% sum

                  }

                })

                R0_hat = R0_hat-used_length
                R0_hat[R0_hat<0] = 0

                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol=N, byrow=T)

                L0_hat = sapply(1:length(R0_hat), function(i) if(sum(unseen_sample[,i]>0)>0) (R0_hat[i] / nrow(unseen)) else 0 )

                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow=T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample)==0)] = 0

                for (i in 1:length(L0_hat)){

                  tip = list(edge=matrix(c(2,1),1,2),
                             tip.label=unseen_name[i],
                             edge.length=L0_hat[i],
                             Nnode=1)
                  class(tip) = "phylo"

                  tree_bt = tree_bt + tip

                }

              } else raw = lapply(raw, function(i) i[1:nrow(data[[1]]),])

            } else {

              p_bt = p_bt[1:nrow(data),]
              raw = lapply(1:ncol(p_bt), function(j){

                lapply(1:nrow(p_bt), function(i) rbinom(n=nT, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)

              })

              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)

            }



            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))

            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)

            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq>0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq>0]

            m_gamma = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_gamma_freq, i, data_type='incidence'))
            m_alpha = sapply(coverage_expected, function(i) coverage_to_size(bootstrap_data_alpha_freq, i, data_type='incidence'))

            aL = phyBranchAL_Inc(phylo = tree_bt, data = bootstrap_data_gamma_raw, datatype = "incidence_raw", rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

            gamma = as.vector(PhD:::PhD.m.est(aL=aL_table_gamma, m=m_gamma, Q=c(0,1,2), datatype='incidence_raw', nt=n) %>% t)


            aL_table_alpha = c()

            for (i in 1:N){

              x = raw[[i]]

              aL = phyBranchAL_Inc(phylo = tree_bt, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

              aL_table_alpha = rbind(aL_table_alpha, aL_table)

            }

            alpha = as.vector((PhD:::PhD.m.est(aL=aL_table_alpha, m=m_alpha, Q=c(0,1,2), datatype='incidence_raw', nt=n)/N) %>% t)

          }

          gamma = gamma#[under_max_alpha]

          alpha = alpha#[under_max_alpha]

          beta = gamma/alpha

          order = rep(c(0,1,2), each=length(coverage_expected))#[under_max_alpha]

          beta = data.frame(Estimate=beta, order)

          C = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(1-order) - 1)/(N^(1-order)-1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(order-1) - 1)/(N^(order-1)-1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate-1)/(N-1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1)))$Estimate

          beta = beta$Estimate

          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along=3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        se = matrix(0, ncol = 7, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)

      }

    }

    if (level=='functional') {

      FD_by_tau = function(data, distance_matrix, tau, coverage_expected, data_type, by) {

        if (data_type=='abundance') {

          zik = data
          zik = zik[rowSums(data)>0,]

          dij = distance_matrix
          dij = dij[rowSums(data)>0, rowSums(data)>0]

          dij[which(dij>tau, arr.ind = T)] = tau
          aik = (1 - dij/tau) %*% as.matrix(zik)
          positive_id = rowSums(aik)>0



          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_v = gamma_x/gamma_a
          # gamma_a = ifelse(gamma_a<1, 1, round(gamma_a))

          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))

          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, c(0,1,2), n) %>% as.vector
          # gamma = FunD:::FD.m.est(ai_vi_gamma, m_gamma, c(0,1,2), n) %>% as.vector


          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          # alpha_a = ifelse(alpha_a<1, 1, round(alpha_a))

          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]

          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))

          if (by=='size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, c(0,1,2), n)/N) %>% as.vector
          if (by=='coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, c(0,1,2), n)/N) %>% as.vector
          # if (by=='size') alpha = (FunD:::FD.m.est(ai_vi_alpha, m_gamma, c(0,1,2), n)/N) %>% as.vector
          # if (by=='coverage') alpha = (FunD:::FD.m.est(ai_vi_alpha, m_alpha, c(0,1,2), n)/N) %>% as.vector

        }

        if (data_type=='incidence_raw') {

          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D

          gamma_Y = data_gamma_freq[-1]

          dij = distance_matrix
          dij = dij[gamma_Y>0, gamma_Y>0]
          gamma_Y = gamma_Y[gamma_Y>0]

          dij[which(dij>tau, arr.ind = T)] = tau
          gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          gamma_a[gamma_a>n] = n

          gamma_v = gamma_Y/gamma_a

          # gamma_a = ifelse(gamma_a<1, 1, round(gamma_a))
          # gamma_a[gamma_a>n] = n

          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))

          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, c(0,1,2), n) %>% as.vector
          # gamma = FunD:::FD.m.est(ai_vi_gamma, m_gamma, c(0,1,2), n) %>% as.vector


          alpha_Y = data_2D[-1,]

          dij = distance_matrix
          dij = dij[rowSums(data_2D[-1,])>0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]

          dij[which(dij>tau, arr.ind = T)] = tau
          alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)

          # alpha_a = ifelse(alpha_a<1, 1, round(alpha_a))
          alpha_a[alpha_a>n] = n
          alpha_a = as.vector(alpha_a)

          # alpha_v = rep(gamma_v, N)
          # alpha_v = alpha_v[alpha_a>0]
          alpha_v = as.vector(as.matrix(alpha_Y))/alpha_a
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]

          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))

          if (by=='size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, c(0,1,2), n)/N) %>% as.vector
          if (by=='coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, c(0,1,2), n)/N) %>% as.vector
          # if (by=='size') alpha = (FunD:::FD.m.est(ai_vi_alpha, m_gamma, c(0,1,2), n)/N) %>% as.vector
          # if (by=='coverage') alpha = (FunD:::FD.m.est(ai_vi_alpha, m_alpha, c(0,1,2), n)/N) %>% as.vector

        }

        return(data.frame(gamma,alpha))

      }

      if (tau_type=='single'){

        if (data_type=='abundance') {

          if (by=='size') {

            output = FD_by_tau(data, distance_matrix, tau, coverage_expected, data_type='abundance', by='size')
            gamma = output$gamma
            alpha = output$alpha

          }

          if (by=='coverage') {

            output = FD_by_tau(data, distance_matrix, tau, coverage_expected, data_type='abundance', by='coverage')
            gamma = output$gamma
            alpha = output$alpha

          }

          gamma = data.frame(coverage_expected, gamma) %>%
            mutate(Method = ifelse(coverage_expected<ref_gamma, 'Interpolated', 'Ref_gamma'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_gamma, m_gamma), 3), Size=rep(m_gamma, 3))

          alpha = data.frame(coverage_expected, alpha) %>%
            mutate(Method = ifelse(coverage_expected<ref_alpha, 'Interpolated',
                                   ifelse(coverage_expected==ref_alpha, 'Ref_alpha',
                                          ifelse(coverage_expected==ref_gamma, 'Ref_gamma', 'Extrapolated'))),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_alpha, m_alpha), 3), Size=rep(m_alpha, 3))

          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha

        }

        if (data_type=='incidence_raw') {

          if (by=='size') {

            output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), distance_matrix, tau, coverage_expected, data_type='incidence_raw', by='size')
            gamma = output$gamma
            alpha = output$alpha

          }

          if (by=='coverage') {

            output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), distance_matrix, tau, coverage_expected, data_type='incidence_raw', by='coverage')
            gamma = output$gamma
            alpha = output$alpha

          }

          gamma = data.frame(coverage_expected, gamma) %>%
            mutate(Method = ifelse(coverage_expected<ref_gamma, 'Interpolated', 'Ref_gamma'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_gamma_freq, m_gamma), 3), Size=rep(m_gamma, 3))

          alpha = data.frame(coverage_expected, alpha) %>%
            mutate(Method = ifelse(coverage_expected<ref_alpha, 'Interpolated',
                                   ifelse(coverage_expected==ref_alpha, 'Ref_alpha',
                                          ifelse(coverage_expected==ref_gamma, 'Ref_gamma', 'Extrapolated'))),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_alpha_freq, m_alpha), 3), Size=rep(m_alpha, 3))

          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha

        }

      }

      if (tau_type=='AUC'){

        cut = seq(0.00000001, 1, length.out = cut_number)
        width = diff(cut)

        if (data_type=='abundance') {

          if (by=='size') {

            gamma_alpha_over_tau = lapply(cut, function(tau) {

              FD_by_tau(data, distance_matrix, tau, coverage_expected, data_type='abundance', by='size')

            })

          }

          if (by=='coverage') {

            gamma_alpha_over_tau = lapply(cut, function(tau) {

              FD_by_tau(data, distance_matrix, tau, coverage_expected, data_type='abundance', by='coverage')

            })

          }
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

          left_limit  = apply(gamma_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

          gamma = colSums((left_limit + right_limit)/2)

          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

          left_limit  = apply(alpha_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

          alpha = colSums((left_limit + right_limit)/2)

          beta_over_tau = gamma_over_tau/alpha_over_tau

          left_limit  = apply(beta_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

          beta = colSums((left_limit + right_limit)/2)

          gamma = data.frame(coverage_expected, gamma) %>%
            mutate(Method = ifelse(coverage_expected<ref_gamma, 'Interpolated', 'Ref_gamma'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_gamma, m_gamma), 3), Size=rep(m_gamma, 3))

          alpha = data.frame(coverage_expected, alpha) %>%
            mutate(Method = ifelse(coverage_expected<ref_alpha, 'Interpolated',
                                   ifelse(coverage_expected==ref_alpha, 'Ref_alpha',
                                          ifelse(coverage_expected==ref_gamma, 'Ref_gamma', 'Extrapolated'))),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_alpha, m_alpha), 3), Size=rep(m_alpha, 3))

          beta = data.frame(coverage_expected, beta) %>%
            mutate(Method = ifelse(coverage_expected<ref_alpha, 'Interpolated',
                                   ifelse(coverage_expected==ref_alpha, 'Ref_alpha',
                                          ifelse(coverage_expected==ref_gamma, 'Ref_gamma', 'Extrapolated'))),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Ind(data_alpha, m_alpha), 3), Size=rep(m_alpha, 3))


        }

        if (data_type=='incidence_raw') {

          if (by=='size') {

            gamma_alpha_over_tau = lapply(cut, function(tau) {

              FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), distance_matrix, tau, coverage_expected, data_type='incidence_raw', by='size')

            })

          }

          if (by=='coverage') {

            gamma_alpha_over_tau = lapply(cut, function(tau) {

              FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), distance_matrix, tau, coverage_expected, data_type='incidence_raw', by='coverage')

            })

          }
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

          left_limit  = apply(gamma_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

          gamma = colSums((left_limit + right_limit)/2)

          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

          left_limit  = apply(alpha_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

          alpha = colSums((left_limit + right_limit)/2)

          beta_over_tau = gamma_over_tau/alpha_over_tau

          left_limit  = apply(beta_over_tau, 1, function(x) x[-cut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

          beta = colSums((left_limit + right_limit)/2)

          gamma = data.frame(coverage_expected, gamma) %>%
            mutate(Method = ifelse(coverage_expected<ref_gamma, 'Interpolated', 'Ref_gamma'),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_gamma_freq, m_gamma), 3), Size=rep(m_gamma, 3))

          alpha = data.frame(coverage_expected, alpha) %>%
            mutate(Method = ifelse(coverage_expected<ref_alpha, 'Interpolated',
                                   ifelse(coverage_expected==ref_alpha, 'Ref_alpha',
                                          ifelse(coverage_expected==ref_gamma, 'Ref_gamma', 'Extrapolated'))),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_alpha_freq, m_alpha), 3), Size=rep(m_alpha, 3))

          beta = data.frame(coverage_expected, beta) %>%
            mutate(Method = ifelse(coverage_expected<ref_alpha, 'Interpolated',
                                   ifelse(coverage_expected==ref_alpha, 'Ref_alpha',
                                          ifelse(coverage_expected==ref_gamma, 'Ref_gamma', 'Extrapolated'))),
                   Order = rep(c(0,1,2), each=length(coverage_expected)/3), Coverage_real=rep(iNEXT:::Chat.Sam(data_alpha_freq, m_alpha), 3), Size=rep(m_alpha, 3))

        }

      }

      gamma = gamma[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      # if (max_alpha_coverage==T) under_max_alpha = !((gamma$Order==0) & (gamma$Coverage_expected>ref_alpha_max)) else under_max_alpha = gamma$Coverage_expected>0
      # gamma = gamma[under_max_alpha,]



      alpha = alpha[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      # alpha = alpha[under_max_alpha,]

      beta = beta[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_expected', 'Coverage_real', 'Size'))

      # beta = beta[under_max_alpha,]

      C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
      U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
      V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
      S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))

      if(nboot>1){

        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"coverage_expected","N",'under_max_alpha',
        #                     'data_type', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))

        # plan(sequential)
        # plan(multiprocess, workers=7)

        # se = parSapply(cl, 1:nboot, function(i){

        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){

          if (data_type=='abundance') {

            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)

            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), distance_matrix, f0_hat, 'abundance')

            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))

            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector

            if (tau_type=='single'){

              if (by=='size') {

                output = FD_by_tau(data_bt, distance_matrix_bt, tau, coverage_expected, data_type='abundance', by='size')
                gamma = output$gamma
                alpha = output$alpha

              }

              if (by=='coverage') {

                output = FD_by_tau(data_bt, distance_matrix_bt, tau, coverage_expected, data_type='abundance', by='coverage')
                gamma = output$gamma
                alpha = output$alpha

              }

              beta=gamma/alpha

            }

            if (tau_type=='AUC'){


              if (by=='size') {

                gamma_alpha_over_tau = lapply(cut, function(tau) {

                  FD_by_tau(data_bt, distance_matrix_bt, tau, coverage_expected, data_type='abundance', by='size')

                })

              }

              if (by=='coverage') {

                gamma_alpha_over_tau = lapply(cut, function(tau) {

                  FD_by_tau(data_bt, distance_matrix_bt, tau, coverage_expected, data_type='abundance', by='coverage')

                })

              }
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              left_limit  = apply(gamma_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

              gamma = colSums((left_limit + right_limit)/2)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              left_limit  = apply(alpha_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

              alpha = colSums((left_limit + right_limit)/2)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              left_limit  = apply(beta_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta = colSums((left_limit + right_limit)/2)

            }

          }

          if (data_type=='incidence_raw') {

            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])

            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), distance_matrix, f0_hat, 'incidence_freq')

            # p_bt = p_bt[rowSums(p_bt)>0,]

            raw = lapply(1:ncol(p_bt), function(j){

              lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)

            })

            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))

            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)

            # data_gamma_freq_bt = data_gamma_freq_bt[data_gamma_freq_bt>0]
            # data_alpha_freq_bt = data_alpha_freq_bt[data_alpha_freq_bt>0]

            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame

            if (tau_type=='single'){

              if (by=='size') {

                output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, coverage_expected, data_type='incidence_raw', by='size')
                gamma = output$gamma
                alpha = output$alpha

              }

              if (by=='coverage') {

                output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, coverage_expected, data_type='incidence_raw', by='coverage')
                gamma = output$gamma
                alpha = output$alpha

              }

              beta = gamma/alpha
            }

            if (tau_type=='AUC'){


              if (by=='size') {

                gamma_alpha_over_tau = lapply(cut, function(tau) {

                  FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, coverage_expected, data_type='incidence_raw', by='size')

                })

              }

              if (by=='coverage') {

                gamma_alpha_over_tau = lapply(cut, function(tau) {

                  FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, coverage_expected, data_type='incidence_raw', by='coverage')

                })

              }
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              left_limit  = apply(gamma_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

              gamma = colSums((left_limit + right_limit)/2)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              left_limit  = apply(alpha_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

              alpha = colSums((left_limit + right_limit)/2)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              left_limit  = apply(beta_over_tau, 1, function(x) x[-cut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta = colSums((left_limit + right_limit)/2)

            }

          }

          # gamma = gamma[under_max_alpha]
          # alpha = alpha[under_max_alpha]
          # beta = beta[under_max_alpha]

          order = rep(c(0,1,2), each=length(coverage_expected))#[under_max_alpha]

          beta = data.frame(Estimate=beta, order)

          C = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(1-order) - 1)/(N^(1-order)-1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(order-1) - 1)/(N^(order-1)-1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate-1)/(N-1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1)))$Estimate

          beta = beta$Estimate

          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along=3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        se = matrix(0, ncol = 7, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)

      }

    }


    se = as.data.frame(se)

    gamma = gamma %>% mutate(LCL = Estimate - tmp * se$gamma,
                             UCL = Estimate + tmp * se$gamma,
                             Region = region_name)

    alpha = alpha %>% mutate(LCL = Estimate - tmp * se$alpha,
                             UCL = Estimate + tmp * se$alpha,
                             Region = region_name)

    beta = beta %>% mutate(  LCL = Estimate - tmp * se$beta,
                             UCL = Estimate + tmp * se$beta,
                             Region = region_name)

    C = C %>% mutate(        LCL = Estimate - tmp * se$C,
                             UCL = Estimate + tmp * se$C,
                             Region = region_name)


    U = U %>% mutate(        LCL = Estimate - tmp * se$U,
                             UCL = Estimate + tmp * se$U,
                             Region = region_name)

    V = V %>% mutate(        LCL = Estimate - tmp * se$V,
                             UCL = Estimate + tmp * se$V,
                             Region = region_name)

    S = S %>% mutate(        LCL = Estimate - tmp * se$S,
                             UCL = Estimate + tmp * se$S,
                             Region = region_name)

    list(gamma = gamma, alpha = alpha, beta = beta, C = C, U = U, V = V, S = S)

  }

  output = lapply(1:length(x), function(i) for_each_region(data = data_list[[i]], region_name = region_names[i], N = Ns[i]))
  names(output) = region_names

  return(output)

}

#' \code{ggiNEXTbeta}: plot the outcome of \code{iNEXTbeta}.
#' @param x the outcome of \code{iNEXTbeta}
#' @param type the plot type: plot of beta diversity (\code{type = "Beta diversity"}) or plot of dissimilarity (\code{type = "Dissimilarity"})
#' @param level the level of diversity computed. only used to decide the label of y axis.
#' @param scale the scale used in \code{facet_grid} in ggplot. Default is \code{'free'}
#' @param main a string describing the title of the plot.
#' @param transp the transparency of the confidence interval. Default is 0.3.
#' @param line_size the size of the lines. Default is 1.
#' @param point_size the size of the points. Default is 1.
#' @param stroke the size of the border of the points. Default is 2.
#' @return a plot of the input outcome of \code{iNEXTbeta}.
#' @import colorRamps
#' @import ggplot2
#' @export
#'
ggiNEXTbeta = function(x, type = c('Beta diversity', 'Dissimilarity'),
                       level = c('Taxonomic', 'Phylogenetic', 'Functional(tau)', 'Functional(AUC)'),
                       scale = 'free', main = NULL, transp = 0.3, line_size = 1, point_size = 1, stroke = 2){

  if (type == 'Beta diversity'){

    gamma = lapply(x, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
    alpha = lapply(x, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
    beta =  lapply(x, function(y) y[["beta"]])  %>% do.call(rbind,.) %>% mutate(div_type = "Beta")  %>% as_tibble()

    df = rbind(gamma, alpha, beta)
    for (i in unique(gamma$Order)) df$Order[df$Order==i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))

    id_obs = which(df$Method == 'Ref_alpha')

    for (i in 1:length(id_obs)) {

      new = df[id_obs[i],]
      new$Coverage_expected = new$Coverage_expected - 0.0001
      new$Method = 'Interpolated'

      newe = df[id_obs[i],]
      newe$Coverage_expected = newe$Coverage_expected + 0.0001
      newe$Method = 'Extrapolated'

      df = rbind(df, new, newe)

    }

    id_obs = which(df$Method == 'Ref_gamma' & df$div_type=='Gamma')

    for (i in 1:length(id_obs)) {

      new = df[id_obs[i],]
      new$Coverage_expected = new$Coverage_expected - 0.0001
      new$Method = 'Interpolated'

      df = rbind(df, new)

    }

    id_obs = which(df$Method == 'Ref_gamma' & df$div_type!='Gamma')

    for (i in 1:length(id_obs)) {

      new = df[id_obs[i],]
      new$Coverage_expected = new$Coverage_expected - 0.0001
      new$Method = 'Extrapolated'

      df = rbind(df, new)

    }

    if (level=='Taxonomic') { ylab = "Taxonomic diversity" }
    if (level=='Phylogenetic') { ylab = "Phylogenetic Hill number" }
    if (level=='Functional(tau)') { ylab = "Functional diversity (given tau)" }
    if (level=='Functional(AUC)') { ylab = "Functional diversity (AUC)" }

  }

  if (type=='Dissimilarity'){

    C = lapply(x, function(y) y[["C"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
    U = lapply(x, function(y) y[["U"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
    V = lapply(x, function(y) y[["V"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-VqN") %>% as_tibble()
    S = lapply(x, function(y) y[["S"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-SqN") %>% as_tibble()

    df = rbind(C, U, V, S)
    for (i in unique(C$Order)) df$Order[df$Order==i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("1-CqN","1-UqN","1-VqN","1-SqN"))

    id_obs = which(df$Method == 'Ref_alpha')

    for (i in 1:length(id_obs)) {

      new = df[id_obs[i],]
      new$Coverage_expected = new$Coverage_expected - 0.0001
      new$Method = 'Interpolated'

      newe = df[id_obs[i],]
      newe$Coverage_expected = newe$Coverage_expected + 0.0001
      newe$Method = 'Extrapolated'

      df = rbind(df, new, newe)

    }

    id_obs = which(df$Method == 'Ref_gamma')

    for (i in 1:length(id_obs)) {

      new = df[id_obs[i],]
      new$Coverage_expected = new$Coverage_expected - 0.0001
      new$Method = 'Extrapolated'

      df = rbind(df, new)

    }

    if (level=='Taxonomic') { ylab = "Taxonomic dissimilarity" }
    if (level=='Phylogenetic') { ylab = "Phylogenetic dissimilarity" }
    if (level=='Functional(tau)') { ylab = "Functional dissimilarity (given tau)" }
    if (level=='Functional(AUC)') { ylab = "Functional dissimilarity (AUC)" }

  }

  lty = c(Interpolated = "solid", Extrapolated = "dotted")
  df$Method = factor(df$Method, levels = c('Interpolated', 'Ref_alpha', 'Extrapolated', 'Ref_gamma'))

  df = mutate(df, region2number = Region)

  for (i in 1:length(unique(df$Region))) df$region2number[df$region2number == unique(df$Region)[i]] = i

  df = mutate(df, point_shape = paste0(Method, '_', region2number))

  sha = c(Ref_alpha_1 = 1, Ref_gamma_1 = 16,
          Ref_alpha_2 = 2, Ref_gamma_2 = 17,
          Ref_alpha_3 = 0, Ref_gamma_3 = 15,
          Ref_alpha_4 = 5, Ref_gamma_4 = 18)

  ggplot(data = df, aes(x = Coverage_expected, y = Estimate, col = Region)) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha=transp) +
    geom_line(data = subset(df, Method %in% c('Interpolated', 'Extrapolated')), aes(linetype=Method), size = line_size) + scale_linetype_manual(values = lty) +
    geom_point(data = subset(df, Method %in% c('Ref_alpha', 'Ref_gamma')), aes(shape = point_shape, stroke = stroke, size = point_size)) +
    scale_shape_manual(values = sha) +
    facet_grid(div_type~Order, scales = scale) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(x='Sample coverage', y=ylab, title=main) +
    guides(shape = FALSE, size=FALSE)
}
