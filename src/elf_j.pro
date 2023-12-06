pro elf_example
  compile_opt idl2
  elf_getspec_crib
  elf_getspec
end

pro write_elf_j_csv, sclet = sclet, zdatatype = zDatatype, trange = trange, suffix = suffix
  compile_opt idl2

  if undefined(zDatatype) then zDatatype = 'pef'
  tstart = trange[0]
  tend = trange[1]
  timeduration = time_double(tend) - time_double(tstart)
  timespan, tstart, timeduration, /seconds ; set the analysis time interval

  ; The low count fluxes can be and should be removed by putting limits on what count level is acceptable based on number of counts.
  elf_init
  elf_load_state, probe = sclet ; Must load ELFIN's position and attitude data the first time
  elf_load_epd, probe = sclet, datatype = zDatatype, type = 'raw'
  elf_getspec, probe = sclet, type = 'raw', /regularize ; now the tplot variables contain raw counts, not eflux

  ; Use these raw counts to produce the df/f error estimate = 1/sqrt(counts) for all quantities you need to. Use calc with globing:
  directions = ['perp', 'omni', 'para', 'anti']
  quant = ['_', zDatatype, '_en_reg_spec2plot_']
  quant = quant.join()

  calc, '"elx' + quant + '????_err"' + '= 1/sqrt(' + '"el' + sclet + quant + '????")' ; <-- what I will use later, err means df/f

  ; reload in eflux units now (two calls: first load data, then compute spectra)
  elf_load_epd, probe = sclet, datatype = zDatatype, type = 'eflux'
  elf_getspec, probe = sclet, type = 'eflux', /regularize ; reload eflux...

  foreach direction, directions do begin
    copy_data, 'el' + sclet + quant + direction, 'elx' + quant + direction
  endforeach

  ; now clean up quantities
  quants2clean = 'elx' + quant + directions ; array of tplot names of energy spectrograms
  errmax2use = 0.5 ; this means % max error is df/f=100*errmax2use
  foreach element, quants2clean do begin
    error2use = element + '_err'
    clean_data, element, error2use, errmax2use
  endforeach

  calc_j_ratios, sclet, quant
  calc_j_ratios, 'x', quant

  directions = [directions, 'paraovrperp']

  if (!true) then begin
    ; NOTE: `quant + directions` gives unreasonable results for some reason
    tplot, 'el' + sclet + quant + directions
    tplot, 'elx' + quant + directions
    filename = 'figures/elx' + quant + suffix

    makepng, filename, /mkdir ; extension appended automatically
  endif

  ; ----------------------------------------------------------------------
  ; Get MLT, L and LAT (IGRF) by tracing to equator
  Re = 6371.0 ; Earth mean radius in km
  tinterpol, 'elx_pos_gsm', 'elx' + quant + 'paraovrperp'
  ttrace2equator, 'elx_pos_gsm_interp', /km
  get_data, 'elx_pos_gsm_interp_foot', data = elx_pos_eq
  L1 = sqrt(total(elx_pos_eq.y ^ 2.0, 2, /nan)) / Re
  store_data, 'elx_L_igrf', data = {x: elx_pos_eq.x, y: L1}

  ; ----------------------------------------------------------------------
  ; Get Local and equotorial magnetic field
  tt89, 'elx_pos_gsm_interp', /igrf_only
  tt89, 'elx_pos_gsm_interp_foot', /igrf_only
  tvectot, 'elx_pos_gsm_interp_bt89', tot = 'elx_Btot'
  tvectot, 'elx_pos_gsm_interp_foot_bt89', tot = 'elx_foot_Btot'
  get_data, 'elx_Btot', data = elx_Btot
  get_data, 'elx_foot_Btot', data = elx_foot_Btot

  ; ----------------------------------------------------------------------
  ; Write state data to csv files
  get_data, 'lossconedeg', data = lossconedeg
  header = ['time', 'alpha_LC', 'L', 'elx_Btot', 'elx_foot_Btot']
  filename = 'data/el' + sclet + '_orb_' + suffix
  write_csv, filename + '.csv', lossconedeg.x, lossconedeg.y, L1, elx_Btot.y, elx_foot_Btot.y, $
    header = header

  ; ----------------------------
  ; Write flux data to csv files
  foreach direction, directions do begin
    get_data, 'elx' + quant + direction, data = elx_en

    filename = 'data/el' + sclet + '_j_' + direction + '_' + suffix
    time = elx_en.x
    j = elx_en.y
    energy = elx_en.v

    id = ncdf_create(filename + '.nc', /clobber)
    xid = ncdf_dimdef(id, 'time', n_elements(time))
    vid = ncdf_dimdef(id, 'energy', n_elements(energy))

    timeid = ncdf_vardef(id, 'time', [xid], /double)
    energyid = ncdf_vardef(id, 'energy', [vid])
    jid = ncdf_vardef(id, 'j', [xid, vid], /double)
    ncdf_control, id, /endef
    ncdf_varput, id, jid, j
    ncdf_varput, id, timeid, time
    ncdf_varput, id, energyid, energy
    ncdf_close, id

    ; write_csv, filename + '.csv', [transpose(time), transpose(j)], $
    ; header = ['time', string(energy)]
  endforeach
end

pro calc_j_ratios, sclet, quant
  ; Calculations related to ratios are refactored to this procedure
  compile_opt idl2
  calc, '"el' + sclet + quant + 'paraovrperp"' + '=' + '"el' + sclet + quant + 'para"' + '/' + '"el' + sclet + quant + 'perp"'
end

pro clean_data, element, error2use, errmax2use
  compile_opt idl2
  ; Refactored code for cleaning quantities to a separate procedure for modularity
  ;
  copy_data, element, 'quant2clean'
  copy_data, error2use, 'error2use'
  get_data, 'quant2clean', data = mydata_quant2clean, dlim = mydlim_quant2clean, lim = mylim_quant2clean
  ntimes = n_elements(mydata_quant2clean.y[*, 0])
  nsectors = n_elements(mydata_quant2clean.y[0, *])
  mydata_quant2clean_temp = reform(mydata_quant2clean.y, ntimes * nsectors)
  get_data, 'error2use', data = mydata_error2use
  mydata_error2use_temp = reform(mydata_error2use.y, ntimes * nsectors)
  ielim = where(abs(mydata_error2use_temp) gt errmax2use, jelim) ; this takes care of NaNs and +/-Inf's as well!
  if jelim gt 0 then mydata_quant2clean_temp[ielim] = !values.f_nan ; make them NaNs, not even zeros
  mydata_quant2clean.y = reform(mydata_quant2clean_temp, ntimes, nsectors) ; back to the original array
  store_data, 'quant2clean', data = mydata_quant2clean, dlim = mydlim_quant2clean, lim = mylim_quant2clean
  copy_data, 'quant2clean', element ; overwrites the previous data in the tplot variable
end

compile_opt idl2
tranges = dictionary( $
  'event01', ['2021-04-17/02:42:00', '2021-04-17/02:46:00'], $
  'event02', ['2021-04-17/04:14:00', '2021-04-17/04:18:00'], $
  'event03', ['2021-04-17/05:47:00', '2021-04-17/05:51:00'] $
  )

sclet = 'b'
foreach trange, tranges, index do begin
  write_elf_j_csv, sclet = sclet, trange = trange, suffix = index
endforeach

end