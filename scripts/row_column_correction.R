row_column_correction<-function(messed_up_dataframe,new_colname_in_qoutes){
	dataset<-messed_up_dataframe
	names(dataset)[-dim(dataset)[2]]<-colnames(dataset)[-1]
	dataset<-dataset[,-dim(dataset)[2]]
	dataset<-cbind(rownames(dataset),dataset)
	names(dataset)[1]<-new_colname_in_qoutes
	return(dataset)	
}

