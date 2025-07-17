
#' Calculate Distance Between Genes and Peaks
#' @description Compute the genomic distances between all genes in the RNA modality and all peaks in the ATAC modality.
#' @param rna A matrix representing the RNA modality data, where columns are features (genes) and rows are cells.
#' @param atac A matrix representing the ATAC modality data, where columns are features (peaks) and rows are cells.
#' @param gene_annotation A data frame containing gene location information (e.g., chromosome, start, end, strand).
#'
#' @return A list containing the modality data and the computed gene-peak distances. This can be used as input for downstream model fitting.
#' @export
gene_peak_distance = function(rna,atac,gene_annotation){

  gene_name = colnames(rna)
  match_gene = match(gene_annotation$gene_name,gene_name)
  gene_loc = gene_annotation[!is.na(match_gene),]  ## Filter the genomic coordinates of genes of interest.
  no_MT = (gene_loc$Chr != "chrM")
  gene_loc = gene_loc[no_MT,]   ## Remove mitochondrial genes

  match_gene_index = match_gene[!is.na(match_gene)]
  index = which(duplicated(match_gene_index)==TRUE)
  match_gene_index = match_gene_index[-index]
  if (length(index) !=0 ) {
     gene_loc = gene_loc[-index,]  #remove duplication
  }

  rna_annot = rna[,gene_loc$gene_name]

  TSS = gene_loc$start        ## Retrieve the positions of transcription start sites (TSS) for genes
  TSS[gene_loc$strand=="-"] = gene_loc$end[gene_loc$strand=="-"]
  gene_location = cbind(gene_loc$Chr,TSS)

  peak_name = colnames(atac)
  peak_name_split = str_split(peak_name,'-',simplify = TRUE)
  peak_location = cbind(peak_name_split[,1],peak_name_split[,2])  ## Retrieve the positions of peaks in the ATAC data


  b_nozero_index = Matrix::Matrix(nrow =dim(peak_location)[1],ncol = dim(gene_location)[1],data=0)
  b_dis = Matrix::Matrix(nrow =dim(peak_location)[1],ncol = dim(gene_location)[1],data=0)
  window_size=200000

  for (i in 1:dim(gene_location)[1]) {
    index_matrix = crossing(i,1:dim(peak_location)[1])
    colnames(index_matrix) = c('gene','peak')
    chr_gene = gene_location[index_matrix$gene,1]
    chr_peak = peak_location[index_matrix$peak,1]
    loc_gene = gene_location[index_matrix$gene,2]
    loc_peak = peak_location[index_matrix$peak,2]
    chr_same = (chr_gene == chr_peak)
    loc_dis = abs(as.numeric(loc_peak) - as.numeric(loc_gene))
    loc_same = abs(as.numeric(loc_peak) - as.numeric(loc_gene)) < window_size
    chr_loc_same = chr_same&loc_same

    b_nozero_index[,i] = as.integer(chr_loc_same)
    b_dis[chr_loc_same,i] = loc_dis[chr_loc_same]
  }


  rna_annot = rna_annot[,Matrix::colSums(b_nozero_index)>0] ## Remove genes without nearby peaks
  atac_annot = atac[,Matrix::rowSums(b_nozero_index)>0]
  b_dis = b_dis[,Matrix::colSums(b_nozero_index)>0]   ## Update the distance matrix
  b_dis = b_dis[Matrix::rowSums(b_nozero_index)>0,]
  b_nozero_index = b_nozero_index[,Matrix::colSums(b_nozero_index)>0]   ## Update the index matrix
  b_nozero_index = b_nozero_index[Matrix::rowSums(b_nozero_index)>0,]
  colnames(b_nozero_index) = colnames(rna_annot)
  rownames(b_nozero_index) = colnames(atac_annot)
  colnames(b_dis) = colnames(rna_annot)
  rownames(b_dis) = colnames(atac_annot)

  return(list(rna_data = rna_annot,atac_data = atac_annot,distance = list(b_01 = b_nozero_index,b_dis = b_dis)))
}


#' Pairwise interactions term
#'
#' @export

pairwise_optim = function(Y, X, B_initial, sigma_init,index_B_01, dist_matrix, lambda_B, rho, w, marginal_par, u,iteration,d0){

  R = (t(Y)%*%Y)/nrow(Y)
  r = (t(X)%*%X)/nrow(X)
  B = B_initial
  sigma_ele = sigma_init


  X_tilde_t = matrix(0,nrow = nrow(r),ncol = ncol(B))
  for (i in 1:ncol(X_tilde_t)) {
    if (sum(index_B_01[,i]==1)==1) {
      X_tilde_t[,i] = matrix(r[,index_B_01[,i]==1],ncol = 1) %*% B[index_B_01[,i]==1,i]
    }else{
      X_tilde_t[,i] = r[,index_B_01[,i]==1] %*% B[index_B_01[,i]==1,i]
    }
  }
  X_tilde<- t(X_tilde_t)

  if (iteration==1) {
    pairwise_iter = 6
  }else{
    pairwise_iter = 5
  }
  for (quadratic_iter in 1:pairwise_iter) {

    for (j in 1:ncol(B)) {
      index = (index_B_01[,j] == 1)

      prior = exp(dist_matrix[index,j]/d0)
      prior = c(prior,0)
      prior = diag(prior**2, nrow = length(prior),ncol = length(prior))

      e = rep(0,nrow(X_tilde))
      e[j] = 1

      if (sum(index)==1) {
        x_tmp = cbind(matrix(X_tilde[, index],ncol = 1),e)
        pairwise_res = pinv(t(x_tmp)%*%x_tmp + (rho/w)*diag(ncol(x_tmp))+ (lambda_B/w)*prior)%*%(t(x_tmp)%*%R[,j]+(rho/w)*(marginal_par[c(index,TRUE),j]+u[c(index,TRUE),j]))

        B[index,j] = pairwise_res[1:(nrow(pairwise_res)-1),]
        sigma_ele[j] = pairwise_res[nrow(pairwise_res),]
        X_tilde[j,]<- t(B[index,j]) %*% matrix(r[index,],nrow=1)

      }else{
        x_tmp = cbind(X_tilde[, index],e)
        pairwise_res = pinv(t(x_tmp)%*%x_tmp + (rho/w)*diag(ncol(x_tmp))+ (lambda_B/w)*prior)%*%(t(x_tmp)%*%R[,j]+(rho/w)*(marginal_par[c(index,TRUE),j]+u[c(index,TRUE),j]))
        B[index,j] = pairwise_res[1:(nrow(pairwise_res)-1),]
        sigma_ele[j] = pairwise_res[nrow(pairwise_res),]
        X_tilde[j,]<- t(B[index,j]) %*% r[index,]
      }
    }

  }
  loss = pairwise_loss(R = R,r = r,B = B,sigma = sigma_ele,dist_matrix = dist_matrix/d0,marginal_par = marginal_par,u = u,lambda_prior = lambda_B,rho = rho,w = w)

  return(list(B,sigma_ele,loss[[1]],loss[[2]]))
}

#' Computing pairwise loss
#'
#' @export
pairwise_loss = function(R,r,B,sigma,dist_matrix,marginal_par,u,lambda_prior,rho,w){
  loss_pairwise = R - diag(sigma) - t(B)%*%r%*%B

  loss_penalty = (lambda_prior)*sum((exp(dist_matrix)*B)**2)
  loss_st = marginal_par - rbind(B,sigma) + u

  loss = w*sum(loss_pairwise^2) + loss_penalty + rho * sum(loss_st^2)
  loss = loss/2
  loss_g = w*sum(loss_pairwise^2)/2 + loss_penalty/2
  return(list(loss,loss_g,loss_penalty/2,rho * sum(loss_st^2)/2,sum((marginal_par - rbind(B,sigma))^2)))
}



#' Nonlinear Least Squares
#'
#' @export
nls_diag = function(Y, X, pairwise_par, marginal_init, u,rho){
  m<- length(Y)
  Y_unpaired_sort<- sort(Y)
  density_Ym = (1:m)/m

  f_xy_function<- function(par,Y_sort_ele,density_Y_ele){
    f_x<- function(y){
      f_x<- 0
      for(i in 1:nrow(X)){
        f_x<- f_x + 1/nrow(X) * 1/(1+exp((X[i,]%*%par[-length(par)]-y)*pi/(3*par[length(par)])^0.5))
      }
      return(f_x)
    }
    f_xy<- ((density_Y_ele - f_x(Y_sort_ele))^2 + (rho/2)/m * sum((par - pairwise_par + u)^2))^0.5
    return(f_xy)
  }

  lower_bound = rep(-Inf,length(pairwise_par))
  lower_bound[length(lower_bound)] = 10^(-10)
  upper_bound = rep(Inf,length(pairwise_par))
  upper_bound[length(upper_bound)] = var(Y)

  model<- nls.lm(par = marginal_init, lower = lower_bound, upper = upper_bound, fn = f_xy_function, Y_sort_ele = Y_unpaired_sort, density_Y_ele = density_Ym,
                 control = nls.lm.control(maxiter=300))
  return(model)
}



#' loss_f
#' @export
loss_f = function(Y, X, parameter){
  m<- length(Y)
  Y_unpaired_sort<- sort(Y)
  density_Ym<- (1:m)/m
  f_x<- function(y){
    f_x<- 0
    for(i in 1:nrow(X)){
      #f_x<- f_x + 1/nrow(X) * 1/(1+exp((X[i,]%*%parameter[-length(parameter)]-y)/parameter[length(parameter)]))
      f_x<- f_x + 1/nrow(X) * 1/(1+exp((X[i,]%*%parameter[-length(parameter)]-y)*pi/(3*parameter[length(parameter)])^0.5))
    }
    return(f_x)
  }
  f_xy<- sum((density_Ym - f_x(Y_unpaired_sort))^2)
  return(f_xy)

}


#' Marginal distribution term
#'
#' @export
marginal_nls = function(Y, X, pairwise_par, marginal_init_des, u, index_B_01_des, nls_diag, loss_f, rho,n){
  foreach(j = 1:n, .combine = rbind, .errorhandling = "pass",.packages = c('bigmemory','minpack.lm'),.export=c('Y', 'X', 'pairwise_par', 'u','index_B_01_des')) %dopar% {
    Y_unpaired <- attach.big.matrix(Y)
    X_unpaired <- attach.big.matrix(X)
    pairwise_parameter <- attach.big.matrix(pairwise_par)
    marginal_init = attach.big.matrix(marginal_init_des)
    u_iter = attach.big.matrix(u)
    index_B_01 <- attach.big.matrix(index_B_01_des)

    if (sum(index_B_01[,j]==1)==1) {
      model_nls <- nls_diag(Y = Y_unpaired[,j],X = matrix(X_unpaired[,index_B_01[,j] == 1],ncol = 1),
                            pairwise_par = pairwise_parameter[c(index_B_01[,j]==1,TRUE),j], marginal_init = marginal_init[c(index_B_01[,j]==1,TRUE),j],
                            u = u_iter[c(index_B_01[,j]==1,TRUE),j],rho = rho)
      loss_value_f = loss_f(Y = Y_unpaired[,j],X = matrix(X_unpaired[,index_B_01[,j] == 1],ncol = 1),parameter = model_nls$par)
    }else{
      model_nls <- nls_diag(Y = Y_unpaired[,j],X = X_unpaired[,index_B_01[,j] == 1],
                            pairwise_par = pairwise_parameter[c(index_B_01[,j]==1,TRUE),j], marginal_init = marginal_init[c(index_B_01[,j]==1,TRUE),j],
                            u = u_iter[c(index_B_01[,j]==1,TRUE),j],rho = rho)
      loss_value_f = loss_f(Y = Y_unpaired[,j],X = X_unpaired[,index_B_01[,j] == 1],parameter = model_nls$par)
    }

    loss_diag = model_nls$deviance
    convergence = model_nls$info
    marginal_par = model_nls$par

    marginal_alpha = rep(0,nrow(index_B_01)+1)
    marginal_alpha[c(index_B_01[,j] == 1,TRUE)] = marginal_par


    c(marginal_alpha,loss_diag,convergence,loss_value_f)
  }
}

#' Main Function of SuperMAP for Learning Cross-Modality Mapping
#'
#' @param input_data The output from the 'gene_peak_distance' function, including both the modality data and the computed gene-peak distances.
#' @param w A hyperparameter that balances the contribution of marginal distributions and pairwise interactions. Default is 2.
#' @param rho A hyperparameter for the ADMM optimization algorithm. Controls the penalty of the constraint term. Default is 1.
#' @param lambda_prior A hyperparameter that controls the confidence level in prior knowledge. Default is 1e-4.
#' @param n_iter Number of iterations to run. Default is 10.
#' @param d0 A hyperparameter defining distance scale, default is 100kb.
#' @param ncore Number of CPU cores to use for parallel computation.
#'
#' @return A list containing the estimated cross-modality mapping and convergence information.
#' @export
supermap = function(input_data, w=2, rho=1, lambda_prior=1*10^(-4), n_iter=10,d0=100000,ncore=60){

  Y_unpaired = input_data$rna_data
  X_unpaired = input_data$atac_data
  index_B_nonzero = as.matrix(input_data$distance$b_01)
  B_dist = as.matrix(input_data$distance$b_dis)

  X_mean = colMeans(X_unpaired)
  Y_mean = colMeans(Y_unpaired)
  X_unpaired = apply(X_unpaired, 2, function(x) x-mean(x))
  Y_unpaired = apply(Y_unpaired, 2, function(x) x-mean(x))


  bigmat_Y_unpaired = as.big.matrix(Y_unpaired)
  bigmat_X_unpaired = as.big.matrix(X_unpaired)
  bigmat_index_B_01 = as.big.matrix(index_B_nonzero)
  descriptor_Y <- describe(bigmat_Y_unpaired)
  descriptor_X <- describe(bigmat_X_unpaired)
  descriptor_index <- describe(bigmat_index_B_01)

  Y_col_variances = apply(Y_unpaired,2,var)
  B_pairwise = index_B_nonzero*0.5
  sigma_pairwise = 0.5*Y_col_variances
  marginal_par = rbind(index_B_nonzero*0.5,0.5*Y_col_variances)
  u_init = matrix(0,nrow = nrow(marginal_par),ncol = ncol(marginal_par))

  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  openblasctl::openblas_set_num_threads(ncore)

  convergence_matrix = matrix(0,nrow = n_iter,ncol = 2)

  for (iteration in 1:n_iter) {
    print(paste0('running iteration ',iteration))
    pairwise_result = pairwise_optim(Y = Y_unpaired,X = X_unpaired,B_initial = B_pairwise,sigma_init = sigma_pairwise,index_B_01 = index_B_nonzero,dist_matrix = B_dist,
                                     lambda_B = lambda_prior,rho = rho,w = w,marginal_par = marginal_par,u = u_init,iteration = iteration,d0 = d0)
    pairwise_par = rbind(pairwise_result[[1]],pairwise_result[[2]])

    bigmat_pairwise_par = as.big.matrix(pairwise_par)
    bigmat_marginal_par = as.big.matrix(marginal_par)
    bigmat_u = as.big.matrix(u_init)
    descriptor_pairwise_par <- describe(bigmat_pairwise_par)
    descriptor_marginal_par <- describe(bigmat_marginal_par)
    descriptor_u <- describe(bigmat_u)
    start_time = Sys.time()
    # marginal_result = marginal_optim(Y = descriptor_Y, X = descriptor_X, pairwise_par = descriptor_pairwise_par,marginal_init_des = descriptor_marginal_par,
    #                                  u = descriptor_u,index_B_01_des = descriptor_index, optim_diag = optim_diag,n=ncol(Y_unpaired),rho = rho,loss_f = loss_f)
    marginal_result = marginal_nls(Y = descriptor_Y, X = descriptor_X, pairwise_par = descriptor_pairwise_par,marginal_init_des = descriptor_marginal_par,
                                   u = descriptor_u,index_B_01_des = descriptor_index, nls_diag = nls_diag,n=ncol(Y_unpaired),rho = rho,loss_f = loss_f)
    print(Sys.time()-start_time)
    marginal_result = t(marginal_result)
    marginal_par = marginal_result[1:nrow(pairwise_par),]
    marginal_loss = marginal_result[nrow(pairwise_par)+1,]
    marginal_convergence = marginal_result[nrow(pairwise_par)+2,]


    u_init = u_init + marginal_par - pairwise_par
    B_pairwise = pairwise_result[[1]]
    sigma_pairwise = pairwise_result[[2]]

    r = sum((marginal_par - pairwise_par)^2)
    g_loss = pairwise_result[[4]]
    f_loss = sum(marginal_result[nrow(pairwise_par)+3,])
    cat(g_loss,f_loss,g_loss+f_loss,r,'\n')
    convergence_matrix[iteration,1] = r
    convergence_matrix[iteration,2] = g_loss+f_loss

    print(paste0('iteration ',iteration,' done'))

  }
  intercept_estimated = as.vector(Y_mean - t(X_mean)%*%B_pairwise)

  stopImplicitCluster()
  stopCluster(cl)
  return(list(estimate_b = B_pairwise,estimated_intercept = intercept_estimated,estimated_sigma = pairwise_result[[2]],
              convergence=convergence_matrix))
}


#' UMAP Visualization of Cells Grouped by Modality
#'
#' @export
cell_type_plot = function(umap,rna_cell_size, rna_label, atac_label){
  rna_umap = as.data.frame(umap[1:rna_cell_size,])
  atac_umap = as.data.frame(umap[(rna_cell_size+1):nrow(umap),])
  rna_umap$cell_type = rna_label[rownames(rna_umap),]
  atac_umap$cell_type = atac_label[rownames(atac_umap),]

  all_cell_types <- unique(rna_umap$cell_type)
  manual_colors <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(all_cell_types)), all_cell_types)

  p = ggplot()+
    geom_point(data=rna_umap,aes(x=umap_1,y = umap_2,color = cell_type),size=0.05)+
    geom_point(data=atac_umap,aes(x=umap_1,y = umap_2,color = cell_type),size=0.5)+
    theme_bw()+
    scale_color_manual(values = manual_colors) +
    labs(x = "UMAP1", y = "UMAP2")+
    theme(panel.grid = element_blank(),#panel.border = element_blank()
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.line = element_line(colour = "black"),legend.position = 'right',
          legend.text = element_text(size = 10),
          axis.title = element_blank(),  # 隐藏坐标轴名称
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1,
          axis.line.x = element_line(color = NA),
          axis.line.y = element_line(color = NA),
          plot.title = element_text(hjust = 0.5))+
    labs(color = "")+
    coord_fixed()
  return(p)
}




#' Performing diagonal integration
#' @export
diagonal_integration = function(supermap_data, learned_mappings = learned_mappings,
                                rna, atac, rna_label){
  b_hat = rbind(learned_mappings$estimated_intercept,learned_mappings$estimate_b)
  atac_data = supermap_data$atac_data
  imputation = cbind(1,atac_data) %*% b_hat

  new_atac = CreateSeuratObject(counts = t(imputation),assay = "ACTIVITY")
  new_atac = NormalizeData(object = new_atac)
  new_atac@assays$ACTIVITY$data = t(imputation)
  new_atac = ScaleData(new_atac)
  gene_use = colnames(imputation)

  transfer.anchors <- FindTransferAnchors(reference = rna, query = new_atac, features = gene_use,scale = FALSE,
                                          reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
  refdata = (rna@assays$RNA$data)[gene_use,]
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,
                             weight.reduction = 'pca',query = new_atac,k.weight = 5)
  cell_label = rna_label
  cell_label = cell_label[colnames(rna),]
  predicted_label = TransferData(anchorset = transfer.anchors, refdata = cell_label,
                                 weight.reduction = 'pca',query = new_atac,k.weight = 5)
  predicted_label = predicted_label@meta.data['predicted.id']
  rna_data = rna@assays$RNA$data[gene_use,]
  atac_data = imputation@assays$id@data[rownames(rna_data),]

  coembed_data <- cbind(rna_data,atac_data)
  coembed = CreateSeuratObject(coembed_data)
  coembed = NormalizeData(object = coembed)
  coembed@assays$RNA$data = coembed_data
  coembed <- ScaleData(coembed, features = gene_use, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = gene_use, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 1:50)
  umap_cord = coembed@reductions$umap@cell.embeddings
  umap_cord = rbind(umap_cord[1:ncol(rna),],deconvol(atac = atac, meta_data = umap_cord))
  pca_cord = coembed@reductions$pca@cell.embeddings
  pca_cord = rbind(pca_cord[1:ncol(rna),],deconvol(atac = atac, meta_data = pca_cord))
  predicted_label = deconvol(atac = atac, meta_data = predicted_label)
  return(list(umap = umap_cord,pca = pca_cord,predicted_label = predicted_label))
}

deconvol = function(atac, meta_data){
  meta_single_correspondence = atac@meta.data['seurat_clusters']
  single_data = matrix(nrow = nrow(meta_single_correspondence),ncol = ncol(meta_data))
  rownames(single_data) = rownames(meta_single_correspondence);colnames(single_data)=colnames(meta_data)
  for (cell in rownames(single_data)) {
    meta_name = meta_single_correspondence[cell,]
    meta_name = paste0("metacell_",meta_name)
    single_data[cell,] = meta_data[meta_name,]
  }
  rownames(single_data) = paste0(rownames(single_data),'_atac')
  return(single_data)
}

#' UMAP Visualization of Cells Grouped by Batch
#' @export
batch_plot = function(umap,rna_cell_size){
  colnames(umap) = c('UMAP_1','UMAP_2')
  p = ggplot()+
    geom_point(data=as.data.frame(umap[1:rna_cell_size,]),aes(x=UMAP_1,y = UMAP_2,color = 'RNA'),size=0.05)+
    geom_point(data=as.data.frame(umap[(rna_cell_size+1):nrow(umap),]),aes(x=UMAP_1,y = UMAP_2,color = 'ATAC'),size=0.5)+
    scale_color_manual(name = "", values = c("RNA" = "#1A77B6", "ATAC" = "#F57F21")) +
    theme_bw()+
    theme(panel.grid = element_blank(),#panel.border = element_blank()
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.line = element_line(colour = "black"),legend.position = 'right',
          legend.text = element_text(size = 10),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1,
          axis.line.x = element_line(color = NA),
          axis.line.y = element_line(color = NA),
          plot.title = element_text(hjust = 0.5))+
    labs(color = "")+
    coord_fixed()
  return(p)
}


#' Visualization of Marker Genes
#' @description Visualize the expression patterns of marker genes across cells.
#' @export
marker_plot = function(cell_umap){
  p = ggplot()+
    geom_point(data=cell_umap,aes_string(x = "umap_1", y = "umap_2", color = "imputed_value"),size=0.1)+
    scale_color_gradient(low = "grey", high = '#CD2626') +
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = 'right',
      legend.text = element_text(size = 13),
      aspect.ratio = 1,
      axis.title = element_text(size = 13)
    )+
    labs(color = "")+
    coord_fixed()+
    guides(color = guide_colorbar(barwidth = 1, barheight = 20))
  return(p)
}

#' Missing Modality Imputation
#'
#' @param supermap_data The output from the 'gene_peak_distance' function, including both the modality data and the computed gene-peak distances.
#' @param learned_mappings Cross-modality mappings learned by the 'supermap' function.
#' @param atac ATAC modality data.
#'
#' @return A matrix containing the imputed values for the missing modality.
#' @export
imputation = function(supermap_data,learned_mappings,atac){
  b_hat = rbind(learned_mappings$estimated_intercept,learned_mappings$estimate_b)
  atac_data = supermap_data$atac_data
  imputation = cbind(1,atac_data) %*% b_hat
  meta_single_correspondence = atac@meta.data['seurat_clusters']
  imputation_single = matrix(nrow = nrow(meta_single_correspondence),ncol = ncol(imputation))
  rownames(imputation_single) = rownames(meta_single_correspondence);colnames(imputation_single)=colnames(imputation)
  for (cell in rownames(imputation_single)) {
    meta_name = meta_single_correspondence[cell,]
    meta_name = paste0("metacell_",meta_name)
    imputation_single[cell,] = imputation[meta_name,]
  }
  return(imputation_single)
}

#' Construct Metacells
#'
#' @param object A Seurat object
#' @param resolution Clustering resolution parameter controlling the granularity of metacell
#' @param graph.name Name of the graph to use for clustering
#' @param algorithm Clustering algorithm to use
#'
#' @return A Seurat object
#' @export
metacell_construct <- function(object, resolution, graph.name, algorithm =1){
  object <- FindClusters(object, graph.name = graph.name,resolution = resolution, algorithm = algorithm)

  # reorder the factor levels
  column <- paste0(graph.name, "_res.", resolution)
  if (!(column %in% colnames(object@meta.data))) column <- "seurat_clusters"
  tmp <- as.character(object@meta.data[[column]])
  levels.clust <- as.character(sort(as.numeric(levels(object@meta.data[[column]]))))
  tmp <- factor(tmp, levels = levels.clust)
  object$seurat_clusters <- tmp
  #if (column != "seurat_clusters") object[[column]] <- NULL
  return(object)
}

#' Compute Metacell-Level Data Matrix
#'
#' @param object A Seurat object
#' @param cluster.name The name of the cluster
#' @param size_factor A scaling factor used to normalize the data
#'
#' @return A matrix representing data at the metacell level.
#' @export
metacell_matrix_RNA <- function(object, cluster.name, size_factor=10^6) {
  clusters <- object@meta.data[[cluster.name]]
  clust.levels <- levels(clusters)
  assay.matrix = object@assays$RNA$counts
  nrow <- length(clust.levels)
  ncol <- nrow(assay.matrix)
  metacell <- matrix(nrow = nrow, ncol = ncol)
  rownames(metacell) <- paste0("metacell_", clust.levels)
  colnames(metacell) <- rownames(assay.matrix)
  for (i in 1:nrow(metacell)) {
    metacell_expression <- apply(assay.matrix[ ,clusters ==as.character(i - 1)], 1, mean)
    metacell[i,] = log(metacell_expression*size_factor/sum(metacell_expression)+1)
  }
  return(metacell)
}

#' Compute Metacell-Level Data Matrix
#'
#' @param object A Seurat object
#' @param cluster.name The name of the cluster
#'
#' @return A matrix representing data at the metacell level.
#' @export
metacell_matrix_ATAC <- function(object, cluster.name) {
  clusters <- object@meta.data[[cluster.name]]
  clust.levels <- levels(clusters)
  assay.matrix = object@assays$ATAC@counts
  assay.matrix = ifelse(as.matrix(assay.matrix) == 0, 0, 1)
  assay.matrix = Matrix(assay.matrix,sparse = TRUE)
  nrow <- length(clust.levels)
  ncol <- nrow(assay.matrix)
  metacell <- matrix(nrow = nrow, ncol = ncol)
  rownames(metacell) <- paste0("metacell_", clust.levels)
  colnames(metacell) <- rownames(assay.matrix)
  for (i in 1:nrow(metacell)) {
    metacell[i, ] <- apply(assay.matrix[, clusters ==as.character(i - 1)], 1, mean)
  }
  return(metacell)
}

#' Compute Metacell-Level Data Matrix
#'
#' @param object A Seurat object
#' @param cluster.name The name of the cluster
#'
#' @return A matrix representing data at the metacell level.
#' @export
metacell_matrix_ADT <- function(object, cluster.name) {
  clusters <- object@meta.data[[cluster.name]]
  clust.levels <- levels(clusters)
  assay.matrix = object@assays$ADT$counts
  nrow <- length(clust.levels)
  ncol <- nrow(assay.matrix)
  metacell <- matrix(nrow = nrow, ncol = ncol)
  rownames(metacell) <- paste0("metacell_", clust.levels)
  colnames(metacell) <- rownames(assay.matrix)
  for (i in 1:nrow(metacell)) {
    metacell_expression <- apply(assay.matrix[ ,clusters ==as.character(i - 1)], 1, mean)
    metacell[i,] = metacell_expression
  }
  return(metacell)
}

