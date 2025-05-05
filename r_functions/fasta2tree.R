fasta2tree <- function(.data_fasta,
                       model_ml = NULL, # new argument: if evolutionary model is provided, skip model test
                       my_outgroup = NULL,
                       n_bootstraps = 100,
                       threshold_bootstraps = 50,
                       n_cpu = parallel::detectCores()){
    
    #### ML Model Selection & Optimization ####
    if(is.null(model_ml)) {
      ## test models, select lowest AIC
      data_modeltest <-
        .data_fasta %>%
        phangorn::modelTest(
          model = "all",
          multicore = ifelse(n_cpu == 1, FALSE, TRUE),
          mc.cores = ifelse(n_cpu == 1, NULL, n_cpu)
        )
      
      #### Initial ModelFit ####
      ## initial fit, use best fit model from prev chunk
      data_modeltest_fit <-
        as.pml(data_modeltest)
      # data_modeltest_fit
      
      # automatically read output from previous line and return model in next line
      modeltest_as.pml_bestfit <- 
        data_modeltest_fit$call$tree %>% 
        as.character() %>% 
        gsub("tree_", "", .)
      
      print(paste("The best evolutionary model is:", modeltest_as.pml_bestfit))
      
    } else {
      
      data_modeltest <-
        .data_fasta %>%
        phangorn::modelTest(
          model = model_ml,
          multicore = ifelse(n_cpu == 1, FALSE, TRUE),
          mc.cores = ifelse(n_cpu == 1, NULL, n_cpu)
        )
      
      #### Initial ModelFit ####
      ## initial fit, use best fit model from prev chunk
      data_modeltest_fit <-
        as.pml(data_modeltest)
      # data_modeltest_fit
      
      # automatically read output from previous line and return model in next line
      modeltest_as.pml_bestfit <- 
        data_modeltest_fit$call$tree %>% 
        as.character() %>% 
        gsub("tree_", "", .)
      
      print(paste("The best version of the ", model_ml, 
                  " evolutionary model is:", modeltest_as.pml_bestfit))
      
    }
    
    #### Optimize Fit ####
    
    ## MAKE SURE TO CHANGE OPTIONS DEPENDING ON YOUR SELECTED MODEL.
    # ?optim.pml gives guidance on how to set optBf and optQ
    # optInv and optGamma should be set depending on whether your model includes +I and/or +G parameters
    
    # Extract the base model by removing +I and +G
    base_model <- 
      gsub("\\+I|\\+G", "", modeltest_as.pml_bestfit) %>%
      trimws()
    
    # Automatically set optBf and optQ based on the base model
    if(base_model == "JC") {
      auto_optBf <- FALSE
      auto_optQ <- FALSE
    } else if(base_model == "K80") {
      auto_optBf <- FALSE
      auto_optQ <- TRUE
    } else if(base_model == "F81") {
      auto_optBf <- TRUE
      auto_optQ <- FALSE
    } else if(base_model == "SYM") {
      auto_optBf <- FALSE
      auto_optQ <- TRUE
    } else {
      # default to optimizing both if the model is not one of the above
      auto_optBf <- TRUE
      auto_optQ <- TRUE
    }
    
    # Determine if the best model includes +I or +G parameters
    auto_optInv <- grepl("\\+I", modeltest_as.pml_bestfit)
    auto_optGamma <- grepl("\\+G", modeltest_as.pml_bestfit)
    
    # print(
    cat("Optimized evolutionary model parameters:\n",
        "optBf =", auto_optBf, "\n",
        "optQ =", auto_optQ, "\n",
        "optInv =", auto_optInv, "\n",
        "optGamma =", auto_optGamma, "\n")
    # )
    
    evolModelFit_opt <-
      data_modeltest_fit %>%
      optim.pml(
        optBf = auto_optBf,
        optQ = auto_optQ,
        optInv = auto_optInv,
        optGamma = auto_optGamma,
        rearrangement = "NNI",
        control = pml.control(trace = 0)
      )
    evolModelFit_opt
    
    #### and Bootstrap ####
    # bootstrap model
    trees_evolModelFit_opt_bs <-
      bootstrap.pml(
        evolModelFit_opt,
        bs = n_bootstraps,
        optNni = TRUE ,
        multicore = ifelse(n_cpu == 1, FALSE, TRUE),
        mc.cores = ifelse(n_cpu == 1, NULL, n_cpu)
      )
    
    ## plotBS functions
    # type = the type of tree to plot, one of "phylogram", "cladogram", "fan", "unrooted", "radial" or "none". If type is "none" the tree is returned with the bootstrap values assigned to the node labels.
    # method = either "FBP" the classical bootstrap (default) or "TBE" (transfer bootstrap)
    # digits = nteger indicating the number of decimal places.
    # p	= only plot support values higher than this percentage number (default is 0).
    
    tree_evolModelFit_opt_bs <-
      plotBS(
        evolModelFit_opt$tree,
        trees_evolModelFit_opt_bs,
        p = threshold_bootstraps/100,
        digits = 2,
        type = "phylogram",
        method = "FBP"
      )
    
    # if max bs val is less than or equal to 1, then convert to scale of 0-100
    conversion_factor <- 
      if (
        max(
          as.numeric(
            tree_evolModelFit_opt_bs$node.label
          ), 
          na.rm = TRUE) <= 1
      ) 100 else 1
    
    # Convert bootstrap support values to values between 0 and 100
    bootstraps_sig <- 
      (as.numeric(tree_evolModelFit_opt_bs$node.label) * conversion_factor) %>% 
      { ifelse(. < threshold_bootstraps, "", .) } %>%
      { ifelse(is.na(.), "", .) }
    
    #### SPECIFY OUTGROUP & ROOT ####
    # this takes the accession number and makes sure that the name of the outgroup matches naming used in tip labels
    
    if(!is.null(my_outgroup)){
      the_outgroup <-
        tree_evolModelFit_opt_bs$tip.label %>%
        as_tibble() %>%
        filter(str_detect(value,
                          gsub("[_\\s].*",
                               "",
                               my_outgroup))) %>%
        pull()
      
      print(
        paste(
          "The outgroup is: ", the_outgroup
        )
      )
    } else {
      the_outgroup <- character()
    }
    
    
    if(length(the_outgroup) > 0) {
      tree_evolModelFit_opt_bs_outgroup <-
        unroot(tree_evolModelFit_opt_bs) %>%
        root(outgroup = the_outgroup)
      
      tree_evolModelFit_opt_bs_outgroup
      
      cat("Outgroup found in the tree. Returning the rerooted tree.\n")
      return(
        list(
          model_ml = data_modeltest_fit,
          best_model = modeltest_as.pml_bestfit,
          model_optim.pml = evolModelFit_opt,
          tree = tree_evolModelFit_opt_bs_outgroup, 
          bootstraps_sig = bootstraps_sig
        )
      )
    } else {
      # If the_outgroup is empty, return the original tree without re-rooting
      cat("Outgroup not found in the tree. Returning the original tree.\n")
      return(
        list(
          model_ml = data_modeltest_fit,
          best_model = modeltest_as.pml_bestfit,
          model_optim.pml = evolModelFit_opt,
          tree = tree_evolModelFit_opt_bs, 
          bootstraps_sig = bootstraps_sig
        )
      )    
    }
  }
