;  $ID: REGRID_GFED.ncl V01 01/12/2014 22:30 BRUCE EXP$
;
;******************************************************************************
;  PROGRAM REGRID_GFED.ncl READS 3-HOURLY EMISSION DATA CALCULATED FROM 
;  GFED_PREPROCESS SUBROUTINE (0.5*0.5), AND REGRID IT INTO WRFCHEM GRID BOX.
;
;  VARIABLES:
;  ============================================================================
;  (1 )
;
;  NOTES:
;  ============================================================================
;  (1 ) ORIGINALLY WRITTEN BY BRUCE. (01/12/2014)
;  (2 ) VARIABLES (iy, im, id) ARE INPUTTED FROM SPECIFICATION OF NCL OPTIONS, 
;       WHICH ARE AS FOLLOWING:
;       ncl iy=2010 im=2 id=1 REGRID_GFED.ncl. (01/13/2014)
;******************************************************************************

;  LOAD FUNCTIONS AND PROCEDURES
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;load "ESMF_regridding.ncl"


begin
; reading the input data (lat-lon grid: GFAS)

; NOTE: this part can be modified according to different cases
    datir_grid     = "../results/grids/"
    datir          = "/Volumes/TOSHIBA_3A/idata/emission/smoke/GFAS/"
    datir_scale3hr = "/Volumes/TOSHIBA_3A/idata/emission/smoke/GFED/GFEDV3/GFED3_3HOURLY/"
    datir_out      = "../results/regrid/GFAS/"

    idom = "d01"   ;"d02"   ; d01 ; domain index
    nx   = 129     ;297     ; 129
    ny   = 84      ;120     ; 84

    mx = 720
    my = 360
    mt = 8

    plot = False ;True
; end modification region
    bv2  = -9999.0

    ; read in the wrf lat and lon info
    rlat = asciiread(datir_grid + "lat_"+idom+".txt", (/ny,nx/), "float")
    rlon = asciiread(datir_grid + "lon_"+idom+".txt", (/ny,nx/), "float")
    ; note: corner_rlat and corner_rlon only for "conserve" interpolation method
    corner_rlat = asciiread(datir_grid + "corner_lat_"+idom+".txt", (/ny,nx,4/), "float")
    corner_rlon = asciiread(datir_grid + "corner_lon_"+idom+".txt", (/ny,nx,4/), "float")

; below to begin the regridding from lat-lon to lambert projection
    srcGridName = "src_SCRIP.nc"
    dstGridName = "dst_SCRIP.nc"
    wgtFileName = "LL_2_CWRF.nc"

    ; for source grids
    Opt                = True
    Opt@LLCorner       = (/  -89.75, -179.75/)
    Opt@URCorner       = (/   89.75,  179.75/)
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    latlon_to_SCRIP(srcGridName,"0.5deg",Opt)
    delete(Opt)

    ; for destination grids
    Opt                = True
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Title          = "WRF Grid"
    Opt@GridCornerLat  = corner_rlat
    Opt@GridCornerLon  = corner_rlon
    curvilinear_to_SCRIP(dstGridName, rlat, rlon, Opt)
    delete(Opt)

    ; for weight
    Opt                      = True
    Opt@InterpMethod         = "conserve"   ; "bilinear","patch","conserve"
    Opt@ForceOverwrite       = True
    Opt@PrintTimings         = True
    Opt@Debug                = True
    ESMF_regrid_gen_weights(srcGridName, dstGridName, wgtFileName, Opt)
    delete(Opt)

    ndays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    isum_days_bk = (/31,59,90,120,151,181,212,243,273,304,334,365/)
    isum_days = new((/12/), "integer")
    isum_days = isum_days_bk 
  
    temp_regrid = new((/mt,ny,nx/), "float")

    dimNames = (/"time", "lat",  "lon"  /)
    dimSizes = (/ -1,       ny,   nx    /)
    dimUnlim = (/True,   False,   False /)

    aernm = (/"OC","BC","PM2_5"/)
    ndate = 365
    aer_in   = new((/3,ndate,my,mx/), "float", bv2)
    gfed_3hr = new((/mt,my,mx/), "float",      bv2)
    aer_3hr  = new((/mt,my,mx/), "float",      bv2)

;   CACULATE JULIAN DAY FROM GREGORIAN DATE FOR THE START DATE OF GFAS DATA
    syear_data  = 2003
    smonth_data = 01
    sday_data   = 01
    sjulian_data = greg2jul(syear_data, smonth_data, sday_data, -1)

;  CALCULATE JULIAN DAY FROM GREGORIAN DATE FOR THE START CACULATING
   sselect_year   = 2009
   sselect_month  = 01
   sselect_day    = 01
   sselect_julian = greg2jul(sselect_year, sselect_month, sselect_day, -1)

;  CALCULATE JULIAN DAY FROM GREGORIAN DATE FOR THE END CACULATING
   eselect_year   = 2009
   eselect_month  = 12
   eselect_day    = 31
   eselect_julian = greg2jul(eselect_year, eselect_month, eselect_day, -1)

; do iy = sselect_year, eselect_year     ;2003,2011

   cyear = iy

   if (iy .eq. 2004 .or. iy .eq. 2008 .or. iy .eq. 2012) then
      ndays(1) = 29
      do in = 1,11
         isum_days(in) = isum_days_bk(in)+1
      end do
   else
      ndays(1) = 28
      do in = 1,11
         isum_days(in) = isum_days_bk(in)
      end do
   end if

;  BASED ON THE START OF GFAS DATA AND THE DATES OF START AND END SELECTION TO
;  CACULATE THE START AND END INDEX OF SELECTION IN ALL THE GFAS DATA
   sindex        = sselect_julian - sjulian_data
   eindex        = eselect_julian - sjulian_data

   do iname = 0,2   ; 0,2 =>"OC", "BC", and "PM2p5"
      src_file = datir + "/GFASv1.0_biomassBurning_" + aernm(iname) + "_2003-2013.nc"
      s        = addfile(src_file, "r")
      aer_in(iname,:,:,:) = s->GFAS(sindex:eindex, :, :)
   end do

;  do im = sselect_month, eselect_month
      if (im .ge. 10) then
         cmonth = im
      else
         cmonth = "0" + im
      end if
    
      src_file1 = datir_scale3hr + cyear + "/fraction_emissions_"+cyear+cmonth+".nc"
      s1        = addfile(src_file1, "r")
      gfed_3hr  = s1->Fraction_of_Emissions

      if (im .eq. 1) then
         julian_day = 0
      else
         julian_day = isum_days(im-2)
      end if

;     do id = 1,ndays(im-1)
         if (id .ge. 10) then
            cday = id
         else
            cday = "0" + id
         end if
   
         julian_day1 = julian_day+id
         print ("juliand_day"+julian_day1)
         do iname = 0,2   ; 0,2 =>"OC", "BC", and "PM2p5"
   
            do ihr = 0,mt-1
               do j = 0,my-1 
                  nj = my-1-j
                  ; unit: g/m2/3hr
                  aer_3hr(ihr,j,:) = aer_in(iname,julian_day1-1,nj,:)*gfed_3hr(ihr,j,:)*1000.*86400.
               end do
            end do
            
            flnm_out = datir_out + cyear + "/GFASv1.0_biomassBurning_"+cyear+cmonth+cday+"_"+aernm(iname)+"_"+idom+".nc"
   
            ; read in the source data (this data from north to south)
            temp = aer_3hr
         
            ; for final data
            Opt                = True
            Opt@PrintTimings   = True
            Opt@ForceOverwrite = True
            Opt@Debug          = True
            do m = 0,mt-1
               temp_regrid(m,:,:) = ESMF_regrid_with_weights(temp(m,:,:),wgtFileName, Opt)
            end do
   
            lat2d = temp_regrid@lat2d
            lon2d = temp_regrid@lon2d
            dim = dimsizes(lat2d)
            nlat = dim(0)
            nlon = dim(1)
   
            system("rm -f "+flnm_out)
   
            n = addfile(flnm_out, "c")
            fAtt               = True            ; assign file attributes
            fAtt@CREATION_DATE = systemfunc ("date")
            fAtt@DECODE_BY     = "ESMF_regridding.ncl was use."
            fileattdef( n, fAtt )            ; copy file attributes
   
            delete(dimNames)
            delete(dimSizes)
            delete(dimUnlim)
            ; define dimension
            dimNames = (/"time", "lat",  "lon"  /)
            dimSizes = (/ -1,       ny,   nx    /)
            dimUnlim = (/True,   False,   False /)
            filedimdef(n,dimNames,dimSizes,dimUnlim)
   
            ; output variables
            filevardef(n, "smoke_emission", "float", (/ "time", "lat", "lon" /) )
            temp_regrid@long_name = "GFAS smoke emission for " + aernm(iname)
            temp_regrid@units     = "g/m2/3hr"
            n->smoke_emission     = temp_regrid
   
            filevardef(n, "lat", "float", (/ "lat", "lon" /) )
            rlat@long_name = "latitude"
            rlat@units     = "degree"
            n->lat        = rlat
   
            filevardef(n, "lon", "float", (/ "lat", "lon" /) )
            rlon@long_name = "longitude"
            rlon@units     = "degree"
            n->lon        = rlon
   
         end do

   if (plot) then
;  ----------------------------------------------------------------------
;   Plot the original and regridded data on a panel plot
;  ----------------------------------------------------------------------
      wks = gsn_open_wks("ps","source")             ; ESMF.ps
      gsn_define_colormap(wks,"gui_default")      ; choose colormap
      plot_orig = new(8,graphic)

      res                     = True              ; Plot mods desired.
      res@gsnDraw             = False             ; We will panel later.
      res@gsnFrame            = False

      res@gsnMaximize         = True              ; Maximize plot
      
      res@mpMaxLonF           = 60  ;max(rlon)
      res@mpMinLonF           = -30 ;min(rlon)
      res@mpMaxLatF           = max(rlat)
      res@mpMinLatF           = min(rlat)

      res@cnFillOn            = True              ; color plot desired
      res@cnLinesOn           = False             ; turn off contour lines
      res@cnLineLabelsOn      = False             ; turn off contour lines

      res@cnLevelSelectionMode=  "ManualLevels"   ; manual levels
      res@cnMinLevelValF      = 0.
      res@cnMaxLevelValF      = 0.0016
      res@cnLevelSpacingF     = 0.0002

      dims = tostring(dimsizes(temp))
      hours = (/"00","03","06","09","12","15","18","21"/)
      do m = 0,mt-1
         res@tiMainString = "GFAS" + hours(m) + "hr smoke emissions: original data (" + \
                            str_join(dims," x ") + ")"
         plot_orig(m)  = gsn_csm_contour_map(wks,temp(m,:,:),res)     
      end do

;  ---Resources for paneling
      pres                  = True
      pres@gsnMaximize      = True
      pres@gsnPanelLabelBar = False

      gsn_panel(wks,plot_orig,(/4,2/),pres)
    
   end if

end 
