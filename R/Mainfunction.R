#' estimates the difference of diversity with two treatments and perform meta analysis
#'
#' \code{iNEXTmeta} is a function that estimates the difference of standardized 3D (taxonomic, phylogenetic and functional) diversity with two treatments (e.g., enhanced vs. control), and perform meta analysis for several studies/sites.
#'
#' @param data (a) For datatype = "abundance", data can be input as a data.frame (species by assemblages). Here an assemblage refers to a combination of study/site and treatment. The names of assemblages must use "_" connect "sites"/"studies" (or other kind of word in your research), site/study and treatment of each assemblage. For example, like "Site_B04_Control" the name of a assemblage, "Site" is the word use in this data, "B04" is the site of this assemblage and "Control" is the treatment of this assemblage.\cr
#' (b) For datatype = "incidence_raw", data can be input as a list with several lists (assemblages) of data.frames, each matrix represents species-by-sampling units incidence data. The names of lists (assemblages) must use "_" connect "sites"/"studies" (or other kind of word in your research), site/study and treatment of each assemblage, as described in datatype "abundance".
#' @param diversity selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, and 'FD' = Functional diversity.
#' @param order.q a number specifying the diversity order, Default is \code{q = 0}.
#' @param datatype data type of input data: individual-based abundance data (datatype = "abundance"), species by sampling-units incidence matrix (datatype = "incidence_raw") with all entries being 0 (non-detection) or 1 (detection).
#' @param base selection of standardization base: sample-size-based (base = "size") or coverage-based (base = "coverage") rarefaction and extrapolation for estimating standardized 3D diversity. Default is \code{base = "coverage"}.
#' @param level a number specifying the particular value of sample coverage (between 0 and 1) or sample sizes that will be used to compute standardized 3D estimates. \cr
#' If base = "coverage" and level = NULL, then this function computes the standardized 3D diversity estimates for the minimum sample coverage among all samples extrapolated to double reference sizes. \cr
#' If base = "size" and level = NULL, then this function computes the standardized 3D diversity estimates for the minimum sample size among all samples extrapolated to double reference sizes.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty for estimating standardized 3D diversity and the associated confidence intervals. Default is 10.
#' @param treatment_order a character vector for the names of treatment. The difference of standardized 3D diversity will be computed as diversity of the first treatment minus the diversity of second treatment.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param PDtree (required only when diversity = "PD"), a phylogenetic tree in Newick format for all observed species in the pooled data.
#' @param PDreftime (required only when diversity = "PD"), a number specifying reference times for PD. Default is NULL (i.e., the age of the root of PDtree).
#' @param PDtype (required only when diversity = "PD"), select PD type: PDtype = "PD" (effective total branch length) or PDtype = "meanPD" (effective number of equally divergent lineages). Default is "meanPD", where meanPD = PD/tree depth.
#' @param FDdistM (required only when diversity = "FD"), a species pairwise distance matrix for all species in the pooled data.
#' @param FDtype (required only when diversity = "FD"), select FD type: FDtype = "tau_values" for FD under specified threshold values, or FDtype = "AUC" (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is "AUC".
#' @param FDtau (required only when diversity = "FD" and FDtype = "tau_values"), a number between 0 and 1 specifying tau values (threshold levels). If NULL (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled data (i.e., quadratic entropy).
#'
#'
#' @import devtools
#' @import iNEXT.3D
#' @import ape
#' @import dplyr
#' @import phyclust
#' @import tidyverse
#' @import stringr
#' @import ade4
#' @import cluster
#'
#'
#' @return a dataframe with columns ’Study/Site (or the word that setting in the names of columns/lists in data)’,
#' ’Order.q’, ’Diversity’, ’Difference’ (difference of diversity for two treatments), ’SE’ (standard deviation of
#' difference), ’LCL’ and ’UCL’. Also, dataframe have two columns with the names of treatments, which are about
#' diversity of two treatments in each study/site. And a column ’weight_fixed’ at the last of dataframe.
#'
#' @examples
#' data("bat")
#' iNEXTmeta(bat$data, diversity="TD", order.q=0,
#'           datatype="incidence_raw", base="coverage", level=0.9,
#'           treatment_order=c("Enhanced","Control"))


iNEXTmeta <- function(data, diversity="TD", order.q=0, datatype="incidence_raw", base="coverage", level=NULL, nboot=10, treatment_order, conf=0.95, PDtree, PDreftime = NULL, PDtype = "meanPD", FDdistM, FDtype = "AUC", FDtau = NULL){


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

  ## each site
  div_T1 <- estimate3D(data_T1, diversity=diversity, q=order.q, datatype=datatype, base=base, nboot=nboot, level=level, PDtree=PDtree, PDreftime=PDreftime, PDtype=PDtype, FDdistM=FDdistM, FDtype=FDtype, FDtau=FDtau)
  div_T2 <- estimate3D(data_T2, diversity=diversity, q=order.q, datatype=datatype, base=base, nboot=nboot, level=level, PDtree=PDtree, PDreftime=PDreftime, PDtype=PDtype, FDdistM=FDdistM, FDtype=FDtype, FDtau=FDtau)

  D1 <- as.data.frame(div_T1)[,5]
  D2 <- as.data.frame(div_T2)[,5]
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
#' \code{ggiNEXTmeta} is a function that provides forest plot for the difference of standardized 3D diversity with two treatments.
#'
#' @param data the outcome of the iNEXTmeta function.
#' @param range the range of the forest plot.
#' @param num_round a number that the values show on the plot are rounded to the specified number of decimal places.
#'
#'
#' @import ggplot2
#' @import GGally
#' @import ggpubr
#' @import forestplot
#'
#'
#' @return a forest plot that visualizing the output of iNEXTmeta.
#'
#'
#' @examples
#' data("bat")
#' output <- iNEXTmeta(bat$data, diversity="TD", order.q=0,
#'                     datatype="incidence_raw", base="coverage", level=0.9,
#'                     treatment_order=c("Enhanced","Control"))
#' ggiNEXTmeta(output, c(-11,16), 3)


ggiNEXTmeta <- function(data,range,num_round){

  name_place <- colnames(data)[1]
  name_treat <- colnames(data)[c(8,9)]
  meta_num <- dim(data)[1]
  name_site <- data[,1][-meta_num]
  diversity <- data$Diversity[1]
  order <- data$Order.q[1]

  forestplot_q <- tibble::tibble(mean  = round(data$Difference[-meta_num],num_round),
                                 lower = round(data$LCL[-meta_num],num_round),
                                 upper = round(data$UCL[-meta_num],num_round),
                                 study = name_site,
                                 q_T1 = round(data[,8][-meta_num],num_round),
                                 q_T2 = round(data[,9][-meta_num],num_round),
                                 diff = round(data$Difference[-meta_num],num_round),
                                 LCL = round(data$LCL[-meta_num],num_round),
                                 UCL= round(data$UCL[-meta_num],num_round),
                                 w_fixed= paste(round(data$weight_fixed[-meta_num],2),"%",sep=""))

  forestplot_q |>
    forestplot(labeltext = c(study, q_T1, q_T2, diff, LCL, UCL, w_fixed),
               clip = range,
               xlog = F) |>
    fp_set_style(box = "royalblue",
                 line = "darkblue",
                 summary = "royalblue",
                 txt_gp = fpTxtGp(ticks = gpar(cex=1))) |>
    fp_add_header(study = c("", name_place),
                  q_T1 = c(paste(diversity,"(q=",order,")",sep=""), name_treat[1]),
                  q_T2 = c(paste(diversity,"(q=",order,")",sep=""), name_treat[2]),
                  diff=c("","Difference"),
                  LCL=c("","LCL"),
                  UCL=c("","UCL"),
                  w_fixed=c("","W(fixed)")) |>
    fp_append_row(mean  = round(data$Difference[meta_num],num_round),
                  lower = round(data$LCL[meta_num],num_round),
                  upper = round(data$UCL[meta_num],num_round),
                  study = "Meta analysis",
                  diff = round(data$Difference[meta_num],num_round),
                  LCL = round(data$LCL[meta_num],num_round),
                  UCL = round(data$UCL[meta_num],num_round),
                  is.summary = TRUE) |>
    fp_decorate_graph(graph.pos = 4) |>
    fp_set_zebra_style("#EFEFEF")
}



