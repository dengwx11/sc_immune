prepare_same_input_in_simulation <- function(same_input, output_directory){
    for(i in seq(length(same_input$X))){
        saveRDS(same_input$X[[i]],paste0(output_directory,"/X_",i,".rds"))
    }
    write.table(same_input$Y0,paste0(output_directory, "/Y0.txt"),quote =F,row.names = F,col.names=F)
    same_input$Y0 <- NULL
    same_input$X <- NULL
    saveRDS(same_input, paste0(output_directory, "/raw_same_input.rds"))
}