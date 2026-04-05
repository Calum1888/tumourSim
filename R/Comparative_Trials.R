km_log_rank_test(control_data, treatment_data){

  test_data <- data.frame(
    time      = c(control_data$time,
                  treatment_data$time),
    status    = c(control_data$status,
                  treatment_data$status),
    treatment = c(rep(0, nrow(control_data)),
                  rep(1, nrow(treatment_data))))

  log_rank_test <- survdiff(Surv(time, status)~treatment,data = test_data)

  p_val <- 1 - pchisq(log_rank_test, df =1)

  return(p_val)

}
