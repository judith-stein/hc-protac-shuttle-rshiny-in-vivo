Model <-rxode2({
  
  # Parameter Variability ---------------------------------------------------------
  
  EC50 <- EC50 * exp(eta.EC50)
  if (EC50 < 0){
    EC50 <- 0
  }
  Ag_cell_t <- Ag_cell_t * exp(eta.Ag_cell_t)
  if (Ag_cell_t < 0){
    Ag_cell_t <- 0
  }

  # Modeling assignemnts -----------------------------------

  SF ~ 10^9 * 1 / 6.023 * 10^-23

  CL_ADC ~ CL_ADC / 24  # convert from days to hours
  CLD_ADC ~ CLD_ADC / 24  # convert from days to hours
  CL_Drug ~ CL_Drug / 24  # convert from days to hours
  CLD_Drug ~ CLD_Drug / 24  # convert from days to hours
  K_ADC_dec ~ K_ADC_dec / 24  # convert from days to hours
  P_ADC ~ P_ADC / 24  # convert from days to hours
  P_Drug ~ P_Drug / 24  # convert from days to hours
  D_ADC ~ D_ADC / 24  # convert from days to hours
  D_Drug ~ D_Drug / 24  # convert from days to hours
  tau ~ tau * 24  # convert from days to hours
  DT_tumor ~ DT_tumor * 24  # convert from days to hours
  CL_Ab ~ CL_Ab / 24  # convert from days to hours
  CLD_Ab ~ CLD_Ab / 24  # convert from days to hours
  P_Ab ~ P_Ab / 24  # convert from days to hours
  D_Ab ~ D_Ab / 24  # convert from days to hours
  CL_Ag ~ CL_Ag / 24  # convert from days to hours
  P_Ag ~ P_Ag / 24  # convert from days to hours
  D_Ag ~ D_Ag / 24  # convert from days to hours

  K_Drug_cyto_dt_on ~ K_Drug_cyto_dt_on * SF / V_cell  # From 1/nM/h to 1/h
  Drug_Target_cell_cyto_t ~ Drug_Target_cell_cyto_t * V_cell / SF  # From nM to 1

  # new approach for multiple dosing for DAR
  DAR <- D / A

  # Add total tumor volume
  V_tumor_mm3 <- V_tumor_pro_mm3 + V_tumor_dyi_1_mm3 + V_tumor_dyi_2_mm3 + V_tumor_dyi_3_mm3
  V_tumor ~ V_tumor_mm3 * 10^-6
  R_Tumor ~ (3 * V_tumor * 1000 / (4 * pi))^(1/3)
  Kg ~ log(2) / DT_tumor * (1 - V_tumor_pro_mm3 / V_tumor_max) / (1 + (log(2) / DT_tumor * V_tumor_pro_mm3 / k_lin)^psi)^(1/psi)

  # Logistic function from Thomas Rysiok
  if (!is.na(Drug_cell_cyto_f) && Drug_cell_cyto_f > 0) {
    nM_Drug_cell_cyto_f ~ Drug_cell_cyto_f * SF / V_cell
  } else {nM_Drug_cell_cyto_f ~ 0}  # OtherBWise log of -1 gets NA
  tl ~ log(nM_Drug_cell_cyto_f) - log(EC50)
  LOGI ~ k_g / (1 + (k_g / k_z - 1) * exp(-k_r * k_g * tl))
  R_Kill ~ k_kill_max * (log(2) / DT_tumor)^f_DT_kill * LOGI
  NC_Tumor ~ NCL_tumor * V_tumor
  if (Kg - R_Kill < 0 && thresholdV_tumor != 0 && V_tumor_pro_mm3 <= thresholdV_tumor) {
    R_Kill ~ 0
    Kg ~ 0
  }


  #####################################

  # Modeling ordinary differential equations --------------------------------

  ## Differential equations
  
  d/dt(Ag_C1_b) <- - (CL_Ag / V_C1_Ag) * Ag_C1_b -
                  (Ag_C1_b / V_C1_Ag - Ag_ex_b / E_Ag) * V_tumor / BW *
                  (2 * P_Ag * R_Cap / R_Krogh^2 + 6 * D_Ag / R_Tumor^2)
                  + K_Ab_cell_ag_on * Ab_C1_f * Ag_C1_f / V_C1_Ab
                  - K_Ab_cell_ag_off * Ag_C1_b
                  + K_ADC_cell_ag_on * ADC_C1_f * Ag_C1_f / V_C1_ADC
                  - K_ADC_cell_ag_off * Ag_C1_b
  
  if (enableTmdd) {
    d/dt(Ag_C1_b) <- d/dt(Ag_C1_b) 
                    + K_Ag_shed * (ADC_TMDD_cell_b_ag + Ab_TMDD_cell_b_ag) * NC_TMDD_cell * SF / BW
  }
  
  d/dt(Ag_C1_f) <- - (CL_Ag / V_C1_Ag) * Ag_C1_f -
                      (Ag_C1_f / V_C1_Ag - Ag_ex_f / E_Ag) * V_tumor / BW *
                      (2 * P_Ag * R_Cap / R_Krogh^2 + 6 * D_Ag / R_Tumor^2)
                    - K_Ab_cell_ag_on * Ab_C1_f * Ag_C1_f / V_C1_Ab
                    + K_Ab_cell_ag_off * Ag_C1_b
                    - K_ADC_cell_ag_on * ADC_C1_f * Ag_C1_f / V_C1_ADC
                    + K_ADC_cell_ag_off * Ag_C1_b
  
  if (enableTmdd) {
    d/dt(Ag_C1_f) <- d/dt(Ag_C1_f) 
                    + K_Ag_shed * (Ag_TMDD_cell_t - ADC_TMDD_cell_b_ag - Ab_TMDD_cell_b_ag) * NC_TMDD_cell * SF / BW
  }
  
  d/dt(Ag_ex_b) <- (Ag_C1_b / V_C1_Ag - Ag_ex_b / E_Ag) *
                    (2 * P_Ag * R_Cap / R_Krogh^2 + 6 * D_Ag / R_Tumor^2)
                  + K_Ab_cell_ag_on * Ab_ex_f / E_Ab * Ag_ex_f / E_Ag
                  - K_Ab_cell_ag_off * Ag_ex_b / E_Ag
                  + K_ADC_cell_ag_on * ADC_ex_f / E_ADC * Ag_ex_f / E_Ag
                  - K_ADC_cell_ag_off * Ag_ex_b / E_Ag
                  + K_Ag_shed * (ADC_cell_b_ag + Ab_cell_b_ag) * NC_Tumor * SF / V_tumor
  
  d/dt(Ag_ex_f) <-  (Ag_C1_f / V_C1_Ag - Ag_ex_f / E_Ag) *
                    (2 * P_Ag * R_Cap / R_Krogh^2 + 6 * D_Ag / R_Tumor^2)
                  - K_Ab_cell_ag_on * Ab_ex_f / E_Ab * Ag_ex_f / E_Ag
                  + K_Ab_cell_ag_off * Ag_ex_b / E_Ag
                  - K_ADC_cell_ag_on * ADC_ex_f / E_ADC * Ag_ex_f / E_Ag
                  + K_ADC_cell_ag_off * Ag_ex_b / E_Ag
                  + K_Ag_shed * (Ag_cell_t - ADC_cell_b_ag - Ab_cell_b_ag) * NC_Tumor * SF / V_tumor
  
  # d/dt(Ag_cell_t) <- - K_ADC_cell_int * (ADC_cell_b_ag + Ab_cell_b_ag)
  #                     - K_Ag_shed * Ag_cell_t
  #                     + K_ADC_cell_int * (ADC_cell_b_ag + Ab_cell_b_ag)
  #                     + K_Ag_shed * Ag_cell_t
  
  # d/dt(Ag_TMDD_cell_t) <- - K_ADC_TMDD_cell_int * (ADC_TMDD_cell_b_ag + Ab_TMDD_cell_b_ag)
  #                     - K_Ag_shed * Ag_TMDD_cell_t
  #                     + K_ADC_TMDD_cell_int * (ADC_TMDD_cell_b_ag + Ab_TMDD_cell_b_ag)
  #                     + K_Ag_shed * Ag_TMDD_cell_t

  d/dt(Ab_C1_f) <- - (CL_Ab / V_C1_Ab) * Ab_C1_f -
    (CLD_Ab / V_C1_Ab) * Ab_C1_f +
    (CLD_Ab / V_C2_Ab) * Ab_C2_f -
    (Ab_C1_f / V_C1_Ab - Ab_ex_f / E_Ab) * V_tumor / BW *
    (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2)
    - K_Ab_cell_ag_on * Ab_C1_f * Ag_C1_f / V_C1_Ab
    + K_Ab_cell_ag_off * Ag_C1_b

  if (enableTmdd) {
    d/dt(Ab_C1_f) <- d/dt(Ab_C1_f) -
      (K_ADC_TMDD_cell_pino * NC_TMDD_cell / (BW * V_C1_Ab) * Ab_C1_f) -
      K_Ab_cell_ag_on * Ab_C1_f * (Ag_TMDD_cell_t - ADC_TMDD_cell_b_ag - Ab_TMDD_cell_b_ag) *
      NC_TMDD_cell * SF / (BW * V_C1_Ab) +
      K_Ab_cell_ag_off * Ab_TMDD_cell_b_ag * NC_TMDD_cell * SF / BW
  }

  d/dt(Ab_C2_f) <- CLD_Ab / V_C1_Ab * Ab_C1_f -
    CLD_Ab / V_C2_Ab * Ab_C2_f

  d/dt(Ab_ex_f) <- (Ab_C1_f / V_C1_Ab - Ab_ex_f / E_Ab) *
    (2 * P_Ab * R_Cap / R_Krogh^2 + 6 * D_Ab / R_Tumor^2) +
    (- K_Ab_cell_ag_on * Ab_ex_f / E_Ab * (Ag_cell_t - ADC_cell_b_ag - Ab_cell_b_ag) +
       K_Ab_cell_ag_off * Ab_cell_b_ag) * NC_Tumor * SF / V_tumor +
    1 / tau * V_tumor_dyi_3_mm3 * 10^5 *
    Ab_cell_b_ag * SF / V_tumor -
    K_ADC_cell_lyso_pino * NCL_tumor / E_Ab * Ab_ex_f
   - K_Ab_cell_ag_on * Ab_ex_f / E_Ab * Ag_ex_f / E_Ag
   + K_Ab_cell_ag_off * Ag_ex_b / E_Ag

  d/dt(Ab_cell_b_ag) <- K_Ab_cell_ag_on * Ab_ex_f / E_Ab * (Ag_cell_t - ADC_cell_b_ag - Ab_cell_b_ag) -
    K_Ab_cell_ag_off * Ab_cell_b_ag - K_ADC_cell_int * Ab_cell_b_ag -
    log(2) / DT_tumor * Ab_cell_b_ag
    - K_Ag_shed * Ab_cell_b_ag

  d/dt(Ab_TMDD_cell_b_ag) <- K_Ab_cell_ag_on * (Ag_TMDD_cell_t - ADC_TMDD_cell_b_ag -
    Ab_TMDD_cell_b_ag) * Ab_C1_f / V_C1_Ab -
    K_Ab_cell_ag_off * Ab_TMDD_cell_b_ag -
    K_ADC_TMDD_cell_int * Ab_TMDD_cell_b_ag
    - K_Ag_shed * Ab_TMDD_cell_b_ag

  d/dt(ADC_C1_f) <- - (CL_ADC / V_C1_ADC) * ADC_C1_f -
    (CLD_ADC / V_C1_ADC) * ADC_C1_f +
    (CLD_ADC / V_C2_ADC) * ADC_C2_f -
    (ADC_C1_f / V_C1_ADC - ADC_ex_f / E_ADC) * V_tumor / BW *
    (2 * P_ADC * R_Cap / R_Krogh^2 + 6 * D_ADC / R_Tumor^2)
    - K_ADC_cell_ag_on * ADC_C1_f * Ag_C1_f / V_C1_ADC
    + K_ADC_cell_ag_off * Ag_C1_b

  if (enableTmdd) {
    d/dt(ADC_C1_f) <- d/dt(ADC_C1_f) -
      (K_ADC_TMDD_cell_pino * NC_TMDD_cell / (BW * V_C1_ADC) * ADC_C1_f) -
      K_ADC_cell_ag_on * ADC_C1_f * (Ag_TMDD_cell_t - ADC_TMDD_cell_b_ag - Ab_TMDD_cell_b_ag) *
      NC_TMDD_cell * SF / (BW * V_C1_ADC) +
      K_ADC_cell_ag_off * ADC_TMDD_cell_b_ag * NC_TMDD_cell * SF / BW
  }

  d/dt(ADC_C2_f) <- CLD_ADC / V_C1_ADC * ADC_C1_f - CLD_ADC / V_C2_ADC * ADC_C2_f

  d/dt(Drug_C1_f) <- - CL_Drug / V_C1_Drug * Drug_C1_f -
    CLD_Drug / V_C1_Drug * Drug_C1_f +
    CLD_Drug / V_C1_Drug * Drug_C2_f +
    (K_ADC_dec * ADC_C1_f * DAR) / V_C1_Drug +
    CL_ADC * DAR * (ADC_C1_f / V_C1_ADC) / V_C1_Drug -
    (Drug_C1_f - Drug_ex_f / (V_tumor * E_Drug)) *
    V_tumor / (V_C1_Drug * BW) *
    (2 * P_Drug * R_Cap / R_Krogh^2 + 6 * D_Drug / R_Tumor^2) -
    K_Drug_ex_ntp_on_off * (1 - f_ex_ub) * Drug_C1_f +
    K_Drug_ex_ntp_on_off * f_ex_ub * Drug_C1_b_ntp

  if (enableTmdd) {
    d/dt(Drug_C1_f) <- d/dt(Drug_C1_f) +
      Drug_TMDD_cell_cyto_f * K_Drug_ex_out * NC_TMDD_cell * SF / (BW * V_C1_Drug) -
      K_Drug_ex_in * V_TMDD_cell * NC_TMDD_cell / (BW * V_C1_Drug) * Drug_C1_f
  }

  d/dt(Drug_C2_f) <- CLD_Drug / V_C2_Drug * Drug_C1_f -
    CLD_Drug / V_C2_Drug * Drug_C2_f

  d/dt(Drug_C1_b_ntp) <- K_Drug_ex_ntp_on_off * (1 - f_ex_ub) * Drug_C1_f - K_Drug_ex_ntp_on_off * f_ex_ub * Drug_C1_b_ntp

  # d/dt(DAR) --> See DAR <- D / A

  d/dt(ADC_ex_f) <- (ADC_C1_f / V_C1_ADC - ADC_ex_f / E_ADC) *
    (2 * P_ADC * R_Cap / R_Krogh^2 + 6 * D_ADC / R_Tumor^2) +
    (- K_ADC_cell_ag_on * ADC_ex_f / E_ADC * (Ag_cell_t - ADC_cell_b_ag - Ab_cell_b_ag) +
       K_ADC_cell_ag_off * ADC_cell_b_ag) * NC_Tumor * SF / V_tumor +
    1 / tau * V_tumor_dyi_3_mm3 * 10^5 *
    (ADC_cell_b_ag + ADC_cell_lyso_f) * SF / V_tumor -
    K_ADC_cell_lyso_pino * NCL_tumor / E_ADC * ADC_ex_f
    - K_ADC_cell_ag_on * ADC_ex_f / E_ADC * Ag_ex_f / E_Ag
    + K_ADC_cell_ag_off * Ag_ex_b / E_Ag

  d/dt(Drug_ex_f) <- (Drug_C1_f - Drug_ex_f / (V_tumor * E_Drug)) * V_tumor *
    (2 * P_Drug * R_Cap / R_Krogh^2 + 6 * D_Drug / R_Tumor^2) +
    K_ADC_dec * ADC_ex_f / E_ADC * DAR * V_tumor +
    (K_ADC_dec * ADC_cell_b_ag * DAR + K_Drug_ex_out * Drug_cell_cyto_f) * NC_Tumor * SF -
    K_Drug_ex_in * NC_Tumor * (V_cell / (V_tumor * E_Drug)) * Drug_ex_f +
    1 / tau * V_tumor_dyi_3_mm3 * 10^5 * (Drug_cell_cyto_f + Drug_cell_cyto_b_dt + Drug_cell_lyso_f) * SF

  d/dt(ADC_cell_b_ag) <- K_ADC_cell_ag_on * ADC_ex_f / E_ADC * (Ag_cell_t - ADC_cell_b_ag - Ab_cell_b_ag) -
    K_ADC_cell_ag_off * ADC_cell_b_ag -
    K_ADC_cell_int * ADC_cell_b_ag -
    log(2) / DT_tumor * ADC_cell_b_ag
    - K_Ag_shed * ADC_cell_b_ag

  d/dt(ADC_cell_lyso_f) <- K_ADC_cell_int * ADC_cell_b_ag -
    K_ADC_deg * ADC_cell_lyso_f +
    K_ADC_cell_lyso_pino * ADC_ex_f / (E_ADC * SF) -
    log(2) / DT_tumor * ADC_cell_lyso_f

  d/dt(Drug_cell_lyso_f) <- K_ADC_deg * ADC_cell_lyso_f * DAR -
    K_Drug_lyso_out * (V_cell / V_cell_lyso) * Drug_cell_lyso_f +
    K_Drug_lyso_in * Drug_cell_cyto_f -
    log(2) / DT_tumor * Drug_cell_lyso_f

  d/dt(Drug_cell_cyto_f) <- K_Drug_lyso_out * (V_cell / V_cell_lyso) * Drug_cell_lyso_f -
    K_Drug_lyso_in * Drug_cell_cyto_f -
    K_Drug_ex_out * Drug_cell_cyto_f -
    K_Drug_cyto_dt_on * Drug_cell_cyto_f * (Drug_Target_cell_cyto_t - Drug_cell_cyto_b_dt) -
    K_Drug_met * Drug_cell_cyto_f +
    K_Drug_cyto_dt_off * Drug_cell_cyto_b_dt +
    K_Drug_ex_in * (V_cell / (V_tumor * E_Drug)) * (Drug_ex_f / SF) -
    log(2) / DT_tumor * Drug_cell_cyto_f

  d/dt(Drug_cell_cyto_b_dt) <- K_Drug_cyto_dt_on * Drug_cell_cyto_f * (Drug_Target_cell_cyto_t - Drug_cell_cyto_b_dt) -
    K_Drug_cyto_dt_off * Drug_cell_cyto_b_dt - log(2) / DT_tumor * Drug_cell_cyto_b_dt

  if (enableTmdd) {
    d/dt(ADC_TMDD_cell_b_ag) <- K_ADC_cell_ag_on * (Ag_TMDD_cell_t - ADC_TMDD_cell_b_ag - Ab_TMDD_cell_b_ag) *
      ADC_C1_f / V_C1_ADC -
      K_ADC_cell_ag_off * ADC_TMDD_cell_b_ag -
      K_ADC_TMDD_cell_int * ADC_TMDD_cell_b_ag
      - K_Ag_shed * ADC_TMDD_cell_b_ag

    d/dt(ADC_TMDD_cell_lyso_f) <- K_ADC_TMDD_cell_int * ADC_TMDD_cell_b_ag -
      K_ADC_TMDD_cell_deg * ADC_TMDD_cell_lyso_f +
      K_ADC_TMDD_cell_pino * ADC_C1_f  / (SF * V_C1_ADC)

    d/dt(Drug_TMDD_cell_lyso_f) <- - K_Drug_lyso_out * (V_TMDD_cell / V_TMDD_cell_lyso) *
      Drug_TMDD_cell_lyso_f +
      K_Drug_lyso_in * Drug_TMDD_cell_cyto_f +
      K_ADC_TMDD_cell_deg * ADC_TMDD_cell_lyso_f * DAR

    d/dt(Drug_TMDD_cell_cyto_f) <- K_Drug_ex_in * V_TMDD_cell * Drug_C1_f / SF -
      K_Drug_ex_out * Drug_TMDD_cell_cyto_f +
      K_Drug_lyso_out * V_TMDD_cell / V_TMDD_cell_lyso * Drug_TMDD_cell_lyso_f -
      K_Drug_lyso_in * Drug_TMDD_cell_cyto_f -
      K_Drug_met * Drug_TMDD_cell_cyto_f
  } else {
    d/dt(ADC_TMDD_cell_b_ag) <- 0
    d/dt(Ab_TMDD_cell_b_ag) <- 0
    d/dt(ADC_TMDD_cell_lyso_f) <- 0
    d/dt(Drug_TMDD_cell_lyso_f) <- 0
    d/dt(Drug_TMDD_cell_cyto_f) <- 0
  }

  d/dt(V_tumor_pro_mm3) <- (Kg - R_Kill) * V_tumor_pro_mm3

  d/dt(V_tumor_dyi_1_mm3) <- R_Kill * V_tumor_pro_mm3 - 1 / tau * V_tumor_dyi_1_mm3

  d/dt(V_tumor_dyi_2_mm3) <- 1 / tau * (V_tumor_dyi_1_mm3 - V_tumor_dyi_2_mm3)

  d/dt(V_tumor_dyi_3_mm3) <- 1 / tau * (V_tumor_dyi_2_mm3 - V_tumor_dyi_3_mm3)

  # neBW approach for multiple dosing for DAR
  d/dt(A) ~ d/dt(ADC_C1_f) + d/dt(ADC_C2_f)
  d/dt(D) ~ - K_ADC_dec * D + D / A * d/dt(A)


  ## Post-definitions:

  d/dt(degADC) <- K_ADC_deg * ADC_cell_lyso_f * NC_Tumor * SF # + ADC_ex_f * V_tumor + ADC_cell_b_ag * NC_Tumor * SF + ADC_cell_lyso_f * NC_Tumor * SF 
})
