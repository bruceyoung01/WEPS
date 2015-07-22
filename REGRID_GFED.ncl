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


begin
; reading the input data (lat-lon grid: GFED= 0.5x0.5)

; NOTE: this part can be modified according to different cases
    datir     = "/Volumes/TOSHIBA_3A/idata/emission/smoke/GFED/GFEDV3/emission_flux_3hourly/"
    datir_grid= "../results/grids/"
    datir_out = "../results/regrid/GFED/"

;  WRFChem GRIDS
    idom = "d01"   ;"d02"   ; d01 ; domain index
    nx   = 129     ;297     ; 129
    ny   = 84      ;120     ; 84

;  GFED GRIDS (GLOBLE)
    mx = 720
    my = 360
    mt = 8

    plot = False     ;True
; end modification region

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

    aernm = (/"OC","BC","PM2p5"/)
    ndays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    temp_regrid = new((/mt,ny,nx/), "float")

    dimNames = (/"time", "lat",  "lon"  /)
    dimSizes = (/ -1,       ny,   nx    /)
    dimUnlim = (/True,   False,   False /)

;do iy = 2010,2010     ;2003,2011
    cyear = iy
;  do im = 1,1           ;1,12
      if (im .ge. 10) then
         cmonth = im
      else
         cmonth = "0" + im
      end if
    
      if (iy .eq. 2004 .or. iy .eq. 2008 .or. iy .eq. 2012) then
         ndays(1) = 29
      else
         ndays(1) = 28
      end if
;     do id = 1,1;ndays(im-1)
         if (id .ge. 10) then
            cday = id
         else
            cday = "0" + id
         end if
   
         do iname = 0,2   ; 0,2 =>"OC", "BC", and "PM2p5"

            flnm = datir + cyear + "/" + aernm(iname) + "/GFED3.1_" + cyear + cmonth + cday + "_" + aernm(iname) + ".txt"
            flnm_out = datir_out + cyear + "/" + aernm(iname) + "/GFED3.1_" + cyear + cmonth + cday + "_" + aernm(iname) + "_" + idom + ".nc"
   
            ; read in the source data (this data from north to south)
            temp = asciiread(flnm, (/mt,my,mx/), "float")
         
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
            temp_regrid@long_name = "GFED smoke emission for " + aernm(iname)
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
;      end do
;   end do
;end do

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
         res@tiMainString = "GFED " + hours(m) + "hr smoke emissions: original data (" + \
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
