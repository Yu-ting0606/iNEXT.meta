#' estimates the difference of diversity with two treatments and perform meta analysis
#'
#' \code{iNEXTmeta} is a function that estimates the difference of standardized 3D (taxonomic, phylogenetic and functional) diversity with two treatments (e.g., enhanced vs. control), and perform meta analysis (fixed effect model) for several studies/sites.
#'
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a data.frame (species by assemblages). Here an assemblage refers to a combination of study/site and treatment. The names of assemblages must use "_" connect "sites"/"studies" (or other kind of word in your research), site/study and treatment of each assemblage. For example, like "Site_B04_Control" the name of a assemblage, "Site" is the word use in this data, "B04" is the site of this assemblage and "Control" is the treatment of this assemblage.\cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a list with several data.frames (assemblages), each data.frames represents species-by-sampling units incidence data. The names of lists (assemblages) must use "_" connect "sites"/"studies" (or other kind of word in your research), site/study and treatment of each assemblage, as described in datatype "abundance".
#' @param diversity selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, and 'FD' = Functional diversity.
#' @param order.q a numerical value specifying the diversity order, Default is \code{q = 0}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param base selection of standardization base: sample-size-based (\code{base = "size"}) or coverage-based (\code{base = "coverage"}) rarefaction and extrapolation for estimating standardized 3D diversity. Default is \code{base = "coverage"}.
#' @param level a numerical value specifying the particular value of sample coverage (between 0 and 1) or sample sizes that will be used to compute standardized 3D estimates. \cr
#' If \code{base = "coverage"} and \code{level = NULL}, then this function computes the standardized 3D diversity estimates for the minimum sample coverage among all samples extrapolated to double reference sizes. \cr
#' If \code{base = "size"} and \code{level = NULL}, then this function computes the standardized 3D diversity estimates for the minimum sample size among all samples extrapolated to double reference sizes.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty for estimating standardized 3D diversity and the associated confidence intervals. Default is 10. If more accurate results are required, set \code{nboot = 100} (or \code{nboot = 200}).
#' @param treatment_order a character vector for the names of treatment. The difference of standardized 3D diversity will be computed as diversity of the first treatment minus the diversity of second treatment.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param PDtree (required only when \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled data.
#' @param PDreftime (required only when \code{diversity = "PD"}), a numerical value specifying reference times for PD. Default is NULL (i.e., the age of the root of PDtree).
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{PDtype = "meanPD"}, where meanPD = PD/tree depth.
#' @param FDdistM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled data.
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under a specified threshold value, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{FDtype = "AUC"}.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical value between 0 and 1 specifying the tau value (threshold level) that will be used to compute FD. If \code{FDtau = "NULL"} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled data (i.e., quadratic entropy).
#' @param FDcut_number (required only when \code{diversity = "FD"} and \code{FDtype = "AUC"}) a numeric number to cut [0, 1] interval into equal-spaced sub-intervals to obtain the AUC value by integrating the tau-profile. Equivalently, the number of tau values that will be considered to compute the integrated AUC value. Default is 30. A larger value can be set to obtain more accurate AUC value.
#'
#'
#'
#' @import iNEXT.3D
#' @import stringr
#'
#'
#' @return a dataframe with columns ’Study/Site (or the word that setting in the names of columns/lists in data)’,
#' ’Order.q’, ’Diversity’(use TD, PD or FD), ’Difference’ (difference of diversity for two treatments), ’SE’ (standard deviation of
#' difference), ’LCL’ and ’UCL’. Also, dataframe have two columns with the names of treatments, which are about
#' diversity of two treatments in each study/site. And a column ’weight_fixed’ (weight of each site/study for fixed effect model) at the last of dataframe.
#'
#' @examples
#'
#' ## Taxonomic diversity for incidence data
#'
#' # Coverage-based standardized TD
#' data("bat_incidence_data")
#' output1c <- iNEXTmeta(data = bat_incidence_data, diversity = "TD", order.q = 0,
#'                       datatype = "incidence_raw", base = "coverage", level = 0.9, nboot = 100,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95)
#' output1c
#'
#' # Sized-based standardized TD
#' data("bat_incidence_data")
#' output1s <- iNEXTmeta(data = bat_incidence_data, diversity = "TD", order.q = 0,
#'                       datatype = "incidence_raw", base = "size", level = 0.9, nboot = 100,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95)
#' output1s
#'
#'
#'
#' ## Phylogenetic diversity for incidence data
#'
#' # Coverage-based standardized PD
#' data("bat_incidence_data")
#' data("bat_tree")
#' output2c <- iNEXTmeta(data = bat_incidence_data, diversity = "PD", order.q = 0,
#'                       datatype = "incidence_raw", base = "coverage", level = 0.9, nboot = 10,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                       PDtree = bat_tree, PDreftime = NULL, PDtype = "meanPD")
#' output2c
#'
#' # Sized-based standardized PD
#' data("bat_incidence_data")
#' data("bat_tree")
#' output2s <- iNEXTmeta(data = bat_incidence_data, diversity = "PD", order.q = 0,
#'                       datatype = "incidence_raw", base = "size", level = 0.9, nboot = 10,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                       PDtree = bat_tree, PDreftime = NULL, PDtype = "meanPD")
#' output2s
#'
#'
#'
#' ## Functional diversity for incidence data
#'
#' # Coverage-based standardized FD
#' data("bat_incidence_data")
#' data("bat_distM")
#' output3c <- iNEXTmeta(data = bat_incidence_data, diversity = "FD", order.q = 0,
#'                       datatype = "incidence_raw", base = "coverage", level = 0.9, nboot = 10,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                       FDdistM = bat_distM, FDtype = "AUC", FDcut_number = 30)
#' output3c
#'
#' # Sized-based standardized FD
#' data("bat_incidence_data")
#' data("bat_distM")
#' output3s <- iNEXTmeta(data = bat_incidence_data, diversity = "FD", order.q = 0,
#'                       datatype = "incidence_raw", base = "size", level = 0.9, nboot = 10,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                       FDdistM = bat_distM, FDtype = "AUC", FDcut_number = 30)
#' output3s
#'
#' @export


iNEXTmeta <- function(data, diversity="TD", order.q=0, datatype="incidence_raw", base="coverage", level=NULL, nboot=10, treatment_order, conf=0.95, PDtree, PDreftime = NULL, PDtype = "meanPD", FDdistM, FDtype = "AUC", FDtau = NULL, FDcut_number=30){


  ##分解名字，取出site,each place和treatment的名稱
  if (datatype=="abundance"){
    site_treat_name <- colnames(data)
  }else if (datatype=="incidence_raw"){
    site_treat_name <- names(data)
  }else{
    stop('Error for datatype.')
  }

  split <- str_split(site_treat_name, pattern = "_", n = Inf, simplify = TRUE)

  if (sum(split=="")!=0){
    stop('The format of column names do not conform to the prescribed form.')
  }

  ### 直接取第一個的地區的用語(batdata用site)，就算別人打錯，每個都不一樣
  name_place <- as.data.frame(split)[1,1]
  name_each_place <- as.data.frame(split)[,2]
  name_treatment <- as.data.frame(split)[,3]
  num_place <- length(name_each_place)

  if (length(unique(name_treatment))!=2){
    stop('There are more than two types of treatment.')
  }
  if (any(sort(unique(name_treatment))!=sort(treatment_order))){
    stop('Treatment_order does not match the column names of data.')
  }
  if (any(table(name_each_place)!=2)){
    stop('There is at least a study/site does not have two treatments or have more than two treatments.')
  }
  compare2 <- sapply(1:length(unique(name_each_place)), function(w) all(sort(unique(name_treatment[name_each_place==unique(name_each_place)[w]]))!=sort(treatment_order)))
  if (any(compare2)){
    stop('There is at least a study/site does not have two treatments that are same as the setting of treatment_order.')
  }

  ##重新排列資料，確保在兩種treatment下的each place順序相同，以便後續直接相減
  if (datatype=="abundance"){

    order_data <- data[,order(factor(name_each_place, levels=unique(name_each_place)))]
    order_name_treatment <- name_treatment[order(factor(name_each_place, levels=unique(name_each_place)))]

    data_T1 <- order_data[,order_name_treatment==treatment_order[1]]
    data_T2 <- order_data[,order_name_treatment==treatment_order[2]]

  }else{

    order_data <- data[order(factor(name_each_place, levels=unique(name_each_place)))]
    order_name_treatment <- name_treatment[order(factor(name_each_place, levels=unique(name_each_place)))]

    data_T1 <- order_data[order_name_treatment==treatment_order[1]]
    data_T2 <- order_data[order_name_treatment==treatment_order[2]]

  }

  CC <- qnorm((1-conf)/2, lower.tail = F)

  ## coverage
  if (datatype == "incidence_raw" & base == "coverage" & is.null(level)){
    Cmax_inc <- c()
    for(i in 1:length(data)){
      Cmax_inc <- c(Cmax_inc, iNEXT.3D:::Coverage(data[[i]], datatype="incidence_raw", 2*ncol(data[[i]])))
    }
    Cmax_inc <- min(Cmax_inc)
    level <- Cmax_inc
  }

  if (datatype == "abundance" & base == "coverage" & is.null(level)){
    Cmax_abun <- c()
    for(i in 1:ncol(data)){
      Cmax_abun <- c(Cmax_abun, iNEXT.3D:::Coverage(data[,i], datatype="abundance", 2*sum(data[,i])))
    }
    Cmax_abun <- min(Cmax_abun)
    level <- Cmax_abun
  }

  if (datatype == "incidence_raw" & base == "size" & is.null(level)){
    Smax_inc <- 2*sapply(1:length(data), function(w) dim(data[[w]])[2])
    Smax_inc <- min(Smax_inc)
    level <- Smax_inc
  }

  if (datatype == "abundance" & base == "size" & is.null(level)){
    Smax_abun <- 2*apply(data,2,sum)
    Smax_abun <- min(Smax_abun)
    level <- Smax_abun
  }

  ## each site
  div_T1 <- estimate3D(data_T1, diversity=diversity, q=order.q, datatype=datatype, base=base, nboot=nboot, level=level, PDtree=PDtree, PDreftime=PDreftime, PDtype=PDtype, FDdistM=FDdistM, FDtype=FDtype, FDtau=FDtau, FDcut_number=30)
  div_T2 <- estimate3D(data_T2, diversity=diversity, q=order.q, datatype=datatype, base=base, nboot=nboot, level=level, PDtree=PDtree, PDreftime=PDreftime, PDtype=PDtype, FDdistM=FDdistM, FDtype=FDtype, FDtau=FDtau, FDcut_number=30)

  D1 <- as.data.frame(div_T1)[,6]
  D2 <- as.data.frame(div_T2)[,6]
  Diff <- D1 - D2
  Var <- div_T1$s.e.^2 + div_T2$s.e.^2
  UCL <- Diff+CC*sqrt(Var)
  LCL <- Diff-CC*sqrt(Var)

  ## fixed model
  den_fixed <- (sum(1/Var))
  delta_fixed <- sum(Diff/Var)/den_fixed


  Totaldata <- data.frame(Site=c(unique(name_each_place),"meta_analysis"),
                          Order.q=rep(order.q,((num_place/2)+1)),
                          Diversity=rep(diversity,((num_place/2)+1)),
                          Difference=c(Diff,delta_fixed),
                          SE=c(sqrt(Var),sqrt(1/den_fixed)),
                          LCL=c(LCL,delta_fixed-CC*sqrt(1/den_fixed)),
                          UCL=c(UCL,delta_fixed+CC*sqrt(1/den_fixed)),
                          T1=as.numeric(c(D1,"")), T2=as.numeric(c(D2,"")),
                          weight_fixed=as.numeric(c(((1/Var)/sum(1/Var))*100,"100")))
  colnames(Totaldata)[c(8,9)] <- treatment_order
  colnames(Totaldata)[1] <- name_place

  return(Totaldata)
}



#' forest plot for the difference of standardized 3D diversity with two treatments
#'
#' \code{ggiNEXTmeta} is a function that provides forest plot for the difference of standardized 3D (taxonomic, phylogenetic and functional) diversity with two treatments.
#'
#' @param output the output of the iNEXTmeta function.
#' @param num_round a numerical value that the values show on the plot are rounded to the specified value of decimal places. Default is 3.
#'
#'
#' @import forestplot
#' @import grid
#'
#'
#' @return a forest plot that visualizing the output of iNEXTmeta. In the plot, it shows the difference of diversity with two treatments for each study/site and meta analysis (fixed effect model).
#'
#'
#' @examples
#' ## Taxonomic diversity for incidence data
#'
#' # Coverage-based standardized TD
#' data("bat_incidence_data")
#' output1c <- iNEXTmeta(data = bat_incidence_data, diversity = "TD", order.q = 0,
#'                       datatype = "incidence_raw", base = "coverage", nboot = 100,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95)
#' ggiNEXTmeta(output1c, num_round = 3)
#'
#' # Sized-based standardized TD
#' data("bat_incidence_data")
#' output1s <- iNEXTmeta(data = bat_incidence_data, diversity = "TD", order.q = 0,
#'                       datatype = "incidence_raw", base = "size", nboot = 100,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95)
#' ggiNEXTmeta(output1s, num_round = 3)
#'
#'
#'
#' ## Phylogenetic diversity for incidence data
#'
#' # Coverage-based standardized PD
#' data("bat_incidence_data")
#' data("bat_tree")
#' output2c <- iNEXTmeta(data = bat_incidence_data, diversity = "PD", order.q = 0,
#'                       datatype = "incidence_raw", base = "coverage", nboot = 10,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                       PDtree = bat_tree, PDreftime = NULL, PDtype = "meanPD")
#' ggiNEXTmeta(output2c, num_round = 3)
#'
#' # Sized-based standardized PD
#' data("bat_incidence_data")
#' data("bat_tree")
#' output2s <- iNEXTmeta(data = bat_incidence_data, diversity = "PD", order.q = 0,
#'                       datatype = "incidence_raw", base = "size", nboot = 10,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                       PDtree = bat_tree, PDreftime = NULL, PDtype = "meanPD")
#' ggiNEXTmeta(output2s, num_round = 3)
#'
#'
#'
#' ## Functional diversity for incidence data
#'
#' # Coverage-based standardized FD
#' data("bat_incidence_data")
#' data("bat_distM")
#' output3c <- iNEXTmeta(data = bat_incidence_data, diversity = "FD", order.q = 0,
#'                       datatype = "incidence_raw", base = "coverage", nboot = 10,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                       FDdistM = bat_distM, FDtype = "AUC", FDcut_number = 30)
#' ggiNEXTmeta(output3c, num_round = 3)
#'
#' # Sized-based standardized FD
#' data("bat_incidence_data")
#' data("bat_distM")
#' output3s <- iNEXTmeta(data = bat_incidence_data, diversity = "FD", order.q = 0,
#'                       datatype = "incidence_raw", base = "size", nboot = 10,
#'                       treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                       FDdistM = bat_distM, FDtype = "AUC", FDcut_number = 30)
#' ggiNEXTmeta(output3s, num_round = 3)
#'
#' @export


ggiNEXTmeta <- function(output, num_round=3){

  name_place <- colnames(output)[1]
  name_treat <- colnames(output)[c(8,9)]
  meta_num <- dim(output)[1]
  name_site <- output[,1][-meta_num]
  diversity <- output$Diversity[1]
  order <- output$Order.q[1]
  range <- c(floor(min(output$LCL)), ceiling(max(output$UCL)))

  forestplot_q <- tibble::tibble(mean  = round(output$Difference[-meta_num],num_round),
                                 lower = round(output$LCL[-meta_num],num_round),
                                 upper = round(output$UCL[-meta_num],num_round),
                                 study = name_site,
                                 q_T1 = round(output[,8][-meta_num],num_round),
                                 q_T2 = round(output[,9][-meta_num],num_round),
                                 diff = round(output$Difference[-meta_num],num_round),
                                 LCL = round(output$LCL[-meta_num],num_round),
                                 UCL= round(output$UCL[-meta_num],num_round),
                                 w_fixed= paste(round(output$weight_fixed[-meta_num],2),"%",sep=""))

  forestplot_q |>
    forestplot(labeltext = c(study, q_T1, q_T2, diff, LCL, UCL, w_fixed),
               clip = range,
               xlog = F, txt_gp = fpTxtGp(cex=1.2, ticks=gpar(cex=1)),
               graphwidth = unit(12, "cm"),
               colgap= unit(3,"mm")) |>
    fp_set_style(box = "royalblue",
                 line = "darkblue",
                 summary = "royalblue") |>
    fp_add_header(study = c("", name_place),
                  q_T1 = c(paste(diversity,"(q=",order,")",sep=""), name_treat[1]),
                  q_T2 = c(paste(diversity,"(q=",order,")",sep=""), name_treat[2]),
                  diff=c("","Difference"),
                  LCL=c("","LCL"),
                  UCL=c("","UCL"),
                  w_fixed=c("","W(fixed)")) |>
    fp_append_row(mean  = round(output$Difference[meta_num],num_round),
                  lower = round(output$LCL[meta_num],num_round),
                  upper = round(output$UCL[meta_num],num_round),
                  study = "Meta analysis",
                  diff = round(output$Difference[meta_num],num_round),
                  LCL = round(output$LCL[meta_num],num_round),
                  UCL = round(output$UCL[meta_num],num_round),
                  is.summary = TRUE) |>
    fp_decorate_graph(graph.pos = 4) |>
    fp_set_zebra_style("#EFEFEF")
}



