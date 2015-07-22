;  $ID: REGRID_SEVIRI.ncl V01 01/12/2014 22:30 BRUCE EXP$
;
;******************************************************************************
;  PROGRAM REGRID_SEVIRI.ncl READS 3-HOURLY EMISSION DATA CALCULATED FROM 
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
;       ncl iy=2010 im=2 id=1 REGRID_SEVIRI.ncl. (01/13/2014)
;******************************************************************************

;  LOAD FUNCTIONS AND PROCEDURES
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"


begin
; reading the input data (lat-lon grid: SEVIRI=1.0x1.0)

; NOTE: this part can be modified according to different cases
    datir     = "/Volumes/TOSHIBA_3A/idata/emission/smoke/FEER-SEVIRI/HOURLY_GRID/"
    datir_out = "/Volumes/TOSHIBA_3A/idata/emission/smoke/FEER-SEVIRI/HOURLY_WRF/test/"
    datir_grid= "../results/grids/"
    datir_out = "../results/regrid/SEVIRI/"

    idom = "d01"   ;"d02"   ; d01 ; domain index
    nx   = 129     ;297     ; 129
    ny   = 84      ;120     ; 84

    mt   = 24
   
    dlon = 1.0
    dlat = 1.0
    mx   = 112
    my   = 82

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

    ndays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    temp_regrid = new((/mt,ny,nx/), "float")

    temp        = new((/my,mx/), "float")
    temp@missingValue = -9999.
    temp@_FillValue   = -9999.

    Opt                = True
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
         flnm_out = datir_out + cyear + "/SEVIRI_Emissions_NSSA_" + cyear + cmonth + cday + "_" + idom + ".nc"
   
         do ihr = 0,mt-1
            if (ihr .ge. 10) then
               chr = ihr
            else
               chr = "0"+ihr
            end if
   
            cdate = cday + chr
            f_log = True

            if (f_log) then
               flnm = datir+cyear+"_data_only/SEVIRI_Emissions_NSSA_"+cyear+cmonth+cdate+"00.dat"

               ;read in the source data (this data from north to south)
               temp = asciiread(flnm, (/my,mx/), "float")

               if (all(ismissing(temp))) then
                temp = 0.0
               end if

               delete(Opt)
               ; for source grids
               Opt                = True
               Opt@LLCorner       = (/  -30.5, -40.5/)
               Opt@URCorner       = (/   50.5,  70.5/)
               Opt@Mask2D         = where(.not.ismissing(temp),1,0)
               Opt@ForceOverwrite = True
               Opt@PrintTimings   = True
               latlon_to_SCRIP(srcGridName,"1deg",Opt)
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

               ; for final data
               Opt                = True
               Opt@PrintTimings   = True
               Opt@ForceOverwrite = True
               Opt@Debug          = True
               temp_regrid(ihr,:,:) = ESMF_regrid_with_weights(temp,wgtFileName, Opt)
            else
               temp_regrid(ihr,:,:) = -9999.
            end if
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
         fAtt@DECODE_BY     = "ESMF_regridding.ncl was used. ALL were made by Dr. Feng Zhang, UNL. zfsu@gmail.com"
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
         temp_regrid@long_name = "SEVIRI Sahel smoke PM2.5 emission"
         temp_regrid@units     = "kg/km2/hr"
         temp_regrid@missingValue = -9999.
         temp_regrid@_FillValue   = -9999.
         n->smoke_emission     = temp_regrid

         filevardef(n, "lat", "float", (/ "lat", "lon" /) )
         rlat@long_name = "latitude"
         rlat@units     = "degree"
         n->lat        = rlat
   
         filevardef(n, "lon", "float", (/ "lat", "lon" /) )
         rlon@long_name = "longitude"
         rlon@units     = "degree"
         n->lon        = rlon
;      end do
;   end do
;end do

end 
