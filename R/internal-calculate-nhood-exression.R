# Monkey patch for miloR::calcNhoodExpression which preserves sparsity, running much faster
.calc_expression <- function(nhoods, data.set, subset.row=NULL, block.size=NULL){
  if (!is.null(subset.row)) {
      data.set <- data.set[subset.row, , drop=FALSE]
  }

  # 1. Normalize neighborhoods while keeping them strictly sparse
  sizes       <- Matrix::colSums(nhoods)
  nhoods.norm <- nhoods %*% Matrix::Diagonal(x = 1 / sizes)
  
  # Coerce to Compressed Sparse Column (dgCMatrix) for optimal C-level operations
  if(!is(nhoods.norm, "dgCMatrix")) {
      nhoods.norm <- as(nhoods.norm, "dgCMatrix")
  }

  # 2. Transpose the sparse neighborhood matrix *once* to get Neighborhoods x Cells
  nhoods.norm_t <- Matrix::t(nhoods.norm)

  # 3. Vectorized Math Trick:
  # Mathematically: Out = Dense %*% Sparse
  # Instead, do: t(Out) = t(Sparse) %*% t(Dense)
  #
  # Shape matching: 
  #   nhoods.norm_t (Neighborhoods x Cells) %*% t(data.set) (Cells x Genes)
  #   = (Neighborhoods x Genes)
  #
  # This triggers the highly optimized `Sparse %*% Dense`
  out_t <- nhoods.norm_t %*% Matrix::t(data.set)

  # 4. Transpose back to original shape (Genes x Neighborhoods)
  out <- Matrix::t(out_t)
  
  rownames(out) <- rownames(data.set)
  colnames(out) <- colnames(nhoods)
  return(out)
}
