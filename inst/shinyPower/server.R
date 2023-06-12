if (exists(".SHINY_MODEL")) {
  mod <- .SHINY_MODEL
} else {
  object <- get(".SHINY_MODEL", envir = unmarked:::.shiny_env)
}

coefs <- unmarked:::check_coefs(NULL, mod, TRUE)

inline_wrap <- function(f, ...){
  out <- f(...)
  div(style='display:inline-block; width: 100px; vertical-align:top', out)
}

get_coef_ui <- function(coefs, nulls=FALSE){

  parbase <- "coef_"
  if(nulls){
    parbase <- "null_"
  }
  out <- list()

  for (i in 1:length(coefs)){
    pars <- coefs[[i]]
    submod_name <- names(coefs)[i]
    inps <- lapply(1:length(pars), function(x){
      par_name <- names(pars)[x]
      inp_name <- paste0(parbase,submod_name,"_",par_name)
      inline_wrap(numericInput, inputId=inp_name, label=par_name,
                  value=0, step=0.01)
    })
    out <- c(out, list(h4(submod_name)), inps)
  }
  out
}

get_coefs <- function(input, nulls=FALSE){
  parbase <- "coef_"
  if(nulls) parbase <- "null_"
  pass <- reactiveValuesToList(input)
  inp_sub <- pass[grepl(parbase,names(pass), fixed=TRUE)]
  inp_sub <- pass[!is.na(names(inp_sub))]
  names(inp_sub) <- gsub(parbase, "", names(inp_sub))
  submods <- gsub("_(.*)$","",names(inp_sub))
  pars <- gsub("^(.*)_","",names(inp_sub))
  out <- lapply(unique(submods), function(x){
    vals <- unlist(inp_sub[which(submods==x)])
    names(vals) <- pars[which(submods==x)]
    vals
  })
  names(out) <- unique(submods)
  out
}

get_design_ui <- function(input, default, name){
  nval <- input[[paste0("ndesign_",name)]]
  inps <- lapply(1:nval, function(x){
      inp_name <- paste0("design_",name,"_",x)
      inline_wrap(numericInput, inputId=inp_name, label=NULL,
                  value=default, min=1, step=1)
    })
  inps
}

get_design <- function(input){
  pass <- reactiveValuesToList(input)
  inp_M <- unlist(pass[grepl("design_sites_",names(pass),fixed=TRUE)])
  inp_M <- inp_M[1:input[["ndesign_sites"]]]
  inp_J <- unlist(pass[grepl("design_obs_",names(pass),fixed=TRUE)])
  inp_J <- inp_J[1:input[["ndesign_obs"]]]
  expand.grid(J=sort(inp_J), M=sort(inp_M), T=1)
  #expand.grid(J=inp_J, M=inp_M, T=1)
}

run_analysis <- function(mod, coefs, alpha, nsim, nulls, design){
  unmarkedPowerList(mod, coefs, design, alpha, nulls, nsim)
}

get_coef_tabset <- function(coefs){
  tabsetPanel(
    tabPanel("Effect sizes", get_coef_ui(coefs)),
    tabPanel("Null hypotheses", get_coef_ui(coefs, nulls=TRUE))
  )
}

get_power_plot <- function(object, param){
  if(inherits(object, "unmarkedPowerList")){
    plot(object, param=param)
  } else {
    plot(1, type="n",xlab="",ylab="",xaxt="n",yaxt="n")
  }
}

get_param_selector <- function(input, object){
  dat <- suppressWarnings(summary(object))
  dat <- dat[dat$M==dat$M[1]&dat$J==dat$J[1]&dat$T==dat$T[1],]
  dat <- dat[dat$Parameter != "(Intercept)",]
  ops <- dat$Parameter
  selectInput("plot_param", "Parameter to plot", choices=ops)
}


function(input, output, session){

  #res_auth <- secure_server(
  #  check_credentials = check_credentials(credentials)
  #)

  #output$auth_output <- renderPrint({
  #  reactiveValuesToList(res_auth)
  #})

  options(unmarked_shiny_session=session)
  output$plot <- renderPlot(plot(mod))
  output$coef_ui <- renderUI(get_coef_tabset(coefs))
  output$coefs <- renderPrint(get_coefs(input))
  output$nulls <- renderPrint(get_coefs(input, nulls=TRUE))
  output$mod <- renderUI(HTML(paste0("<b>Model:</b> ","mod")))
  output$class <- renderUI(HTML(paste0("<b>Type:</b>&nbsp&nbsp&nbsp",
                                       class(mod)[1])))
  output$sites <- renderUI(HTML(paste0("<b>Sites:</b>&nbsp&nbsp",
                                       numSites(mod@data))))
  output$design_sites <- renderUI(get_design_ui(input,numSites(mod@data),"sites"))
  output$design_obs <- renderUI(get_design_ui(input,obsNum(mod@data),"obs"))

  observeEvent(input$run, {
    coefs <- isolate(get_coefs(input))
    nulls <- isolate(get_coefs(input, nulls=TRUE))
    design <- isolate(get_design(input))
    alpha <- isolate(input$alpha)
    nsims <- isolate(input$nsims)
    pa <- run_analysis(mod, coefs, alpha, nsims, nulls, design)
    output$summary <- renderTable(
                      suppressWarnings(summary(pa))
                      )
    output$param_selector <- renderUI(get_param_selector(input, pa))
    output$plot <- renderPlot(suppressWarnings(get_power_plot(pa, input$plot_param)))
  })
}



