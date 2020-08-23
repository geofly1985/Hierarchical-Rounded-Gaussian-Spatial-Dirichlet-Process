
cal_stat= function(glcm){

  res=matrix(0,nrow=nrow(glcm),ncol=5)

  for (i in 1:nrow(glcm)) {

  tmp=glcm[i,]/sum(glcm[i,])
  ###########  features ###########
  r=rep(1:sqrt(ncol(glcm)),times=sqrt(ncol(glcm)))
  c=rep(1:sqrt(ncol(glcm)),each=sqrt(ncol(glcm)))

  uX = sum(r*tmp)
  uY = sum(c*tmp)

  sX = sum((r-uX)^2*tmp)
  sY = sum((c-uY)^2*tmp)


  contrast= sum(abs(r - c)^2*tmp)

  correlation= (sum(r*c*tmp)- uX*uY)/sqrt(sX*sY)

  energy= sum(tmp^2)

  entropy= -sum(tmp*log(tmp+10^(-10)))

  homogeneity= sum(tmp/( 1 + (r - c)^2))

  res[i,]=c(contrast,correlation,energy,entropy,homogeneity)
  }

  colnames(res)=c('contrast','correlation','energy','entropy','homogeneity')

  res
}




error_rate=function(mat){

  max=apply(mat,1,max)

  if (ncol(mat) > 5){

    col_pos=c(which(mat[1,]==max[1]),which(mat[2,]==max[2]),which(mat[3,]==max[3]),which(mat[4,]==max[4]),which(mat[5,]==max[5]))

    a11=mat[which(mat==max(mat),arr.ind=TRUE)[1,][1],which(mat==max(mat),arr.ind=TRUE)[1,][2]]

    mat1=mat[-which(mat==max(mat),arr.ind=TRUE)[1,][1],-which(mat==max(mat),arr.ind=TRUE)[1,][2]]
    a22=mat1[which(mat1==max(mat1),arr.ind=TRUE)[1,][1],which(mat1==max(mat1),arr.ind=TRUE)[1,][2]]

    mat2=mat1[-which(mat1==max(mat1),arr.ind=TRUE)[1,][1],-which(mat1==max(mat1),arr.ind=TRUE)[1,][2]]
    a33=mat2[which(mat2==max(mat2),arr.ind=TRUE)[1,][1],which(mat2==max(mat2),arr.ind=TRUE)[1,][2]]

    mat3=mat2[-which(mat2==max(mat2),arr.ind=TRUE)[1,][1],-which(mat2==max(mat2),arr.ind=TRUE)[1,][2]]
    a44=mat3[which(mat3==max(mat3),arr.ind=TRUE)[1,][1],which(mat3==max(mat3),arr.ind=TRUE)[1,][2]]

    mat4 = mat3[-which(mat3==max(mat3),arr.ind=TRUE)[1,][1],-which(mat3==max(mat3),arr.ind=TRUE)[1,][2]]
    a55 = max(mat4)

    err=1.0-(a11+a22+a33+a44+a55)/100

  } else if (ncol(mat) == 5){
    col_pos=c(which(mat[1,]==max[1]),which(mat[2,]==max[2]),which(mat[3,]==max[3]),which(mat[4,]==max[4]),which(mat[5,]==max[5]))

    true_cluster=unique(col_pos)

    if (length(true_cluster)==5){

      err=1.0-(mat[1,col_pos[1]]+mat[2,col_pos[2]]+mat[3,col_pos[3]]+mat[4,col_pos[4]]+mat[5,col_pos[5]])/100

    } else {

      a11=mat[which(mat==max(mat),arr.ind=TRUE)[1,][1],which(mat==max(mat),arr.ind=TRUE)[1,][2]]

      mat1=mat[-which(mat==max(mat),arr.ind=TRUE)[1,][1],-which(mat==max(mat),arr.ind=TRUE)[1,][2]]
      a22=mat1[which(mat1==max(mat1),arr.ind=TRUE)[1,][1],which(mat1==max(mat1),arr.ind=TRUE)[1,][2]]

      mat2=mat1[-which(mat1==max(mat1),arr.ind=TRUE)[1,][1],-which(mat1==max(mat1),arr.ind=TRUE)[1,][2]]
      a33=mat2[which(mat2==max(mat2),arr.ind=TRUE)[1,][1],which(mat2==max(mat2),arr.ind=TRUE)[1,][2]]

      mat3=mat2[-which(mat2==max(mat2),arr.ind=TRUE)[1,][1],-which(mat2==max(mat2),arr.ind=TRUE)[1,][2]]
      a44=mat3[which(mat3==max(mat3),arr.ind=TRUE)[1,][1],which(mat3==max(mat3),arr.ind=TRUE)[1,][2]]

      a55=mat3[-which(mat3==max(mat3),arr.ind=TRUE)[1,][1],-which(mat3==max(mat3),arr.ind=TRUE)[1,][2]]

      err=1.0-(a11+a22+a33+a44+a55)/100

    }
  } else if (ncol(mat) == 4){

    col_pos=c(which(mat[1,]==max[1]),which(mat[2,]==max[2]),which(mat[3,]==max[3]),which(mat[4,]==max[4]),which(mat[5,]==max[5]))

    a11=mat[which(mat==max(mat),arr.ind=TRUE)[1,][1],which(mat==max(mat),arr.ind=TRUE)[1,][2]]

    mat1=mat[-which(mat==max(mat),arr.ind=TRUE)[1,][1],-which(mat==max(mat),arr.ind=TRUE)[1,][2]]
    a22=mat1[which(mat1==max(mat1),arr.ind=TRUE)[1,][1],which(mat1==max(mat1),arr.ind=TRUE)[1,][2]]

    mat2=mat1[-which(mat1==max(mat1),arr.ind=TRUE)[1,][1],-which(mat1==max(mat1),arr.ind=TRUE)[1,][2]]
    a33=mat2[which(mat2==max(mat2),arr.ind=TRUE)[1,][1],which(mat2==max(mat2),arr.ind=TRUE)[1,][2]]

    mat3=mat2[-which(mat2==max(mat2),arr.ind=TRUE)[1,][1],-which(mat2==max(mat2),arr.ind=TRUE)[1,][2]]
    a44= max(mat3)

    err=1.0-(a11+a22+a33+a44)/100

  } else if (ncol(mat) == 3){

    col_pos=c(which(mat[1,]==max[1]),which(mat[2,]==max[2]),which(mat[3,]==max[3]),which(mat[4,]==max[4]),which(mat[5,]==max[5]))

    a11=mat[which(mat==max(mat),arr.ind=TRUE)[1,][1],which(mat==max(mat),arr.ind=TRUE)[1,][2]]

    mat1=mat[-which(mat==max(mat),arr.ind=TRUE)[1,][1],-which(mat==max(mat),arr.ind=TRUE)[1,][2]]
    a22=mat1[which(mat1==max(mat1),arr.ind=TRUE)[1,][1],which(mat1==max(mat1),arr.ind=TRUE)[1,][2]]

    mat2=mat1[-which(mat1==max(mat1),arr.ind=TRUE)[1,][1],-which(mat1==max(mat1),arr.ind=TRUE)[1,][2]]
    a33=max(mat2)

    err=1.0-(a11+a22+a33)/100


  } else {

    col_pos=c(which(mat[1,]==max[1]),which(mat[2,]==max[2]),which(mat[3,]==max[3]),which(mat[4,]==max[4]),which(mat[5,]==max[5]))

    a11=mat[which(mat==max(mat),arr.ind=TRUE)[1,][1],which(mat==max(mat),arr.ind=TRUE)[1,][2]]

    mat1=mat[-which(mat==max(mat),arr.ind=TRUE)[1,][1],-which(mat==max(mat),arr.ind=TRUE)[1,][2]]
    a22=max(mat1)

    err=1.0-(a11+a22)/100


  }

  err
}
