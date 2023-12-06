pro test_v2, ratio = ratio
  compile_opt idl2
  elf_init
  sclet = 'b'

  ; Orbit1
  tstart = '2021-04-17/02:42:00'
  tend = '2021-04-17/02:44:00'
  time2plot = [tstart, tend]
  timeduration = time_double(tend) - time_double(tstart)
  timespan, tstart, timeduration, /seconds ; set the analysis time interval

  elf_load_state, probe = sclet ; you must load ELFIN's position and attitude data the first time (probe='a' which is also the default)
  elf_load_epd, probe = sclet, datatype = 'pef' ; you must first load the EPD data, here, for further use (default is ELFIN-A 'a', 'e', and 'nflux'
  elf_getspec, /regularize, probe = sclet, datatype = 'pef' ; get some spectra (default is 'a',species='e' (or datatype='pef'),type='nflux'  <- note these MUST match the previous call

  calc, ' "el?_pef_en_reg_spec2plot_paraovrperp" = "el?_pef_en_reg_spec2plot_para" / "el?_pef_en_reg_spec2plot_perp" '
  calc, ' "el?_pef_en_reg_spec2plot_antiovrperp" = "el?_pef_en_reg_spec2plot_anti" / "el?_pef_en_reg_spec2plot_perp" '

  if (keyword_set(ratio)) then begin
    get_data, 'elb_pef_en_reg_spec2plot_paraovrperp', data = d
    flux_ratio = mean(d.y, dimension = 1, /nan)
    flux_ratio_min = min(d.y, dimension = 1, /nan)
    flux_ratio_max = max(d.y, dimension = 1, /nan)
    energy_range = d.v
    myPlot = plot(energy_range, flux_ratio, color = 'Blue')
    myPlot.xlog = 1
    myPlot.ylog = 1
    myPlot.title = 'ELFIN-B EPDE precipitating/trapped flux ratio, 2021-04-17'
    myPlot.xtitle = 'Energy [keV]'
    myPlot.ytitle = 'Flux_ratio'
  endif else begin
    get_data, 'elb_pef_en_reg_spec2plot_para', data = d
    flux = mean(d.y, dimension = 1, /nan)
    energy_range = d.v
    myPlot = plot(energy_range, flux, color = 'Blue')
    myPlot.xlog = 1
    myPlot.ylog = 1
    myPlot.title = 'ELFIN-B EPDE precipitating flux, 2021-04-17'
    myPlot.xtitle = 'Energy [keV]'
    myPlot.ytitle = 'Flux'
  endelse
  stop

  ; Orbit2
  tstart = '2021-04-17/04:15:00'
  tend = '2021-04-17/04:16:00'
  time2plot = [tstart, tend]
  timeduration = time_double(tend) - time_double(tstart)
  timespan, tstart, timeduration, /seconds ; set the analysis time interval

  elf_load_state, probe = sclet ; you must load ELFIN's position and attitude data the first time (probe='a' which is also the default)
  elf_load_epd, probe = sclet, datatype = 'pef' ; you must first load the EPD data, here, for further use (default is ELFIN-A 'a', 'e', and 'nflux'
  elf_getspec, /regularize, probe = sclet, datatype = 'pef' ; get some spectra (default is 'a',species='e' (or datatype='pef'),type='nflux'  <- note these MUST match the previous call

  calc, ' "el?_pef_en_reg_spec2plot_paraovrperp" = "el?_pef_en_reg_spec2plot_para" / "el?_pef_en_reg_spec2plot_perp" '
  calc, ' "el?_pef_en_reg_spec2plot_antiovrperp" = "el?_pef_en_reg_spec2plot_anti" / "el?_pef_en_reg_spec2plot_perp" '

  if (keyword_set(ratio)) then begin
    get_data, 'elb_pef_en_reg_spec2plot_paraovrperp', data = d
    flux_ratio = mean(d.y, dimension = 1, /nan)
    energy_range = d.v
    myPlot = plot(energy_range, flux_ratio, color = 'Green', linestyle = ':', /overplot)
  endif else begin
    get_data, 'elb_pef_en_reg_spec2plot_para', data = d
    flux = mean(d.y, dimension = 1, /nan)
    energy_range = d.v
    myPlot = plot(energy_range, flux, color = 'Green', linestyle = ':', /overplot)
  endelse
  stop

  ; Orbit3
  tstart = '2021-04-17/05:48:00'
  tend = '2021-04-17/05:49:00'
  time2plot = [tstart, tend]
  timeduration = time_double(tend) - time_double(tstart)
  timespan, tstart, timeduration, /seconds ; set the analysis time interval

  elf_load_state, probe = sclet ; you must load ELFIN's position and attitude data the first time (probe='a' which is also the default)
  elf_load_epd, probe = sclet, datatype = 'pef' ; you must first load the EPD data, here, for further use (default is ELFIN-A 'a', 'e', and 'nflux'
  elf_getspec, /regularize, probe = sclet, datatype = 'pef' ; get some spectra (default is 'a',species='e' (or datatype='pef'),type='nflux'  <- note these MUST match the previous call

  calc, ' "el?_pef_en_reg_spec2plot_paraovrperp" = "el?_pef_en_reg_spec2plot_para" / "el?_pef_en_reg_spec2plot_perp" '
  calc, ' "el?_pef_en_reg_spec2plot_antiovrperp" = "el?_pef_en_reg_spec2plot_anti" / "el?_pef_en_reg_spec2plot_perp" '

  if (keyword_set(ratio)) then begin
    get_data, 'elb_pef_en_reg_spec2plot_paraovrperp', data = d
    flux_ratio = mean(d.y, dimension = 1, /nan)
    energy_range = d.v
    myPlot = plot(energy_range, flux_ratio, color = 'Red', linestyle = '--', /overplot)
  endif else begin
    get_data, 'elb_pef_en_reg_spec2plot_para', data = d
    flux = mean(d.y, dimension = 1, /nan)
    energy_range = d.v
    myPlot = plot(energy_range, flux, color = 'Red', linestyle = '--', /overplot)
  endelse
  stop

  ; Orbit4 (Elfin A)
  sclet = 'a'
  tstart = '2021-04-17/05:16:00'
  tend = '2021-04-17/05:17:00'
  time2plot = [tstart, tend]
  timeduration = time_double(tend) - time_double(tstart)
  timespan, tstart, timeduration, /seconds ; set the analysis time interval

  elf_load_state, probe = sclet ; you must load ELFIN's position and attitude data the first time (probe='a' which is also the default)
  elf_load_epd, probe = sclet, datatype = 'pef' ; you must first load the EPD data, here, for further use (default is ELFIN-A 'a', 'e', and 'nflux'
  elf_getspec, /regularize, probe = sclet, datatype = 'pef' ; get some spectra (default is 'a',species='e' (or datatype='pef'),type='nflux'  <- note these MUST match the previous call

  calc, ' "el?_pef_en_reg_spec2plot_paraovrperp" = "el?_pef_en_reg_spec2plot_para" / "el?_pef_en_reg_spec2plot_perp" '
  calc, ' "el?_pef_en_reg_spec2plot_antiovrperp" = "el?_pef_en_reg_spec2plot_anti" / "el?_pef_en_reg_spec2plot_perp" '

  if (keyword_set(ratio)) then begin
    get_data, 'ela_pef_en_reg_spec2plot_paraovrperp', data = d
    flux_ratio = mean(d.y, dimension = 1, /nan)
    energy_range = d.v
    myPlot = plot(energy_range, flux_ratio, color = 'Cyan', linestyle = '-', /overplot)
  endif else begin
    get_data, 'ela_pef_en_reg_spec2plot_para', data = d
    flux = mean(d.y, dimension = 1, /nan)
    energy_range = d.v
    myPlot = plot(energy_range, flux, color = 'Cyan', linestyle = '-', /overplot)
  endelse
  stop

  if (keyword_set(ratio)) then begin
    text01 = text(0.2, 0.80, 'Orbit #01', color = !color.blue, /current)
    text02 = text(0.2, 0.75, 'Orbit #02', color = !color.green, /current)
    text03 = text(0.2, 0.70, 'Orbit #03', color = !color.red, /current)
    text04 = text(0.2, 0.65, 'Orbit #04', color = !color.cyan, /current)
    myPlot.Save, 'flux_ratio.png'
  endif else begin
    text01 = text(0.2, 0.80, 'Orbit #01', color = !color.blue, /current)
    text02 = text(0.2, 0.75, 'Orbit #02', color = !color.green, /current)
    text03 = text(0.2, 0.70, 'Orbit #03', color = !color.red, /current)
    text04 = text(0.2, 0.65, 'Orbit #04', color = !color.cyan, /current)
    myPlot.Save, 'flux.png'
  endelse
  stop
end

pro test_v1
  compile_opt idl2
  elf_init
  sclet = 'b'

  tstart = '2021-04-17/02:42:00'
  tend = '2021-04-17/02:44:00'
  filename = 'elb_20210417_024200'

  time2plot = [tstart, tend]
  timeduration = time_double(tend) - time_double(tstart)
  timespan, tstart, timeduration, /seconds ; set the analysis time interval

  elf_load_state, probe = sclet ; you must load ELFIN's position and attitude data the first time (probe='a' which is also the default)
  elf_load_epd, probe = sclet, datatype = 'pef' ; you must first load the EPD data, here, for further use (default is ELFIN-A 'a', 'e', and 'nflux'
  elf_getspec, /regularize, probe = sclet, datatype = 'pef' ; get some spectra (default is 'a',species='e' (or datatype='pef'),type='nflux'  <- note these MUST match the previous call
  ;

  calc, ' "el?_pef_en_reg_spec2plot_paraovrperp" = "el?_pef_en_reg_spec2plot_para" / "el?_pef_en_reg_spec2plot_perp" '
  calc, ' "el?_pef_en_reg_spec2plot_antiovrperp" = "el?_pef_en_reg_spec2plot_anti" / "el?_pef_en_reg_spec2plot_perp" '

  zlim, 'el?_pef_en_reg_spec2plot_paraovrperp', 0.001, 1., 1.
  zlim, 'el?_pef_en_reg_spec2plot_antiovrperp', 0.001, 1., 1.

  tplot, '*en_reg_spec2plot*'
  makepng, filename

  stop

  tstart = '2021-04-17/04:15:00'
  tend = '2021-04-17/04:16:00'
  filename = 'elb_20210417_041500'

  time2plot = [tstart, tend]
  timeduration = time_double(tend) - time_double(tstart)
  timespan, tstart, timeduration, /seconds ; set the analysis time interval

  elf_load_state, probe = sclet ; you must load ELFIN's position and attitude data the first time (probe='a' which is also the default)
  elf_load_epd, probe = sclet, datatype = 'pef' ; you must first load the EPD data, here, for further use (default is ELFIN-A 'a', 'e', and 'nflux'
  elf_getspec, /regularize, probe = sclet, datatype = 'pef' ; get some spectra (default is 'a',species='e' (or datatype='pef'),type='nflux'  <- note these MUST match the previous call
  ;

  calc, ' "el?_pef_en_reg_spec2plot_paraovrperp" = "el?_pef_en_reg_spec2plot_para" / "el?_pef_en_reg_spec2plot_perp" '
  calc, ' "el?_pef_en_reg_spec2plot_antiovrperp" = "el?_pef_en_reg_spec2plot_anti" / "el?_pef_en_reg_spec2plot_perp" '

  zlim, 'el?_pef_en_reg_spec2plot_paraovrperp', 0.001, 1., 1.
  zlim, 'el?_pef_en_reg_spec2plot_antiovrperp', 0.001, 1., 1.

  tplot, '*en_reg_spec2plot*'
  makepng, filename

  stop

  tstart = '2021-04-17/05:48:00'
  tend = '2021-04-17/05:49:00'
  filename = 'elb_20210417_054800'

  time2plot = [tstart, tend]
  timeduration = time_double(tend) - time_double(tstart)
  timespan, tstart, timeduration, /seconds ; set the analysis time interval

  elf_load_state, probe = sclet ; you must load ELFIN's position and attitude data the first time (probe='a' which is also the default)
  elf_load_epd, probe = sclet, datatype = 'pef' ; you must first load the EPD data, here, for further use (default is ELFIN-A 'a', 'e', and 'nflux'
  elf_getspec, /regularize, probe = sclet, datatype = 'pef' ; get some spectra (default is 'a',species='e' (or datatype='pef'),type='nflux'  <- note these MUST match the previous call
  ;

  calc, ' "el?_pef_en_reg_spec2plot_paraovrperp" = "el?_pef_en_reg_spec2plot_para" / "el?_pef_en_reg_spec2plot_perp" '
  calc, ' "el?_pef_en_reg_spec2plot_antiovrperp" = "el?_pef_en_reg_spec2plot_anti" / "el?_pef_en_reg_spec2plot_perp" '

  zlim, 'el?_pef_en_reg_spec2plot_paraovrperp', 0.001, 1., 1.
  zlim, 'el?_pef_en_reg_spec2plot_antiovrperp', 0.001, 1., 1.

  tplot, '*en_reg_spec2plot*'
  makepng, filename

  stop
end

compile_opt idl2
test_v2
; test_v2, /ratio
end