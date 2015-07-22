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
load "/Users/bruce/Documents/A/program/ncl/ncl_lib/function/int2string_format.ncl"


begin
;   reading the input data (lat-lon grid: QFED=0.1x0.1 AND 0.3125x0.25)
;   NOTE: this part can be modified according to different cases
    datir_scale3hr = "/Volumes/TOSHIBA_3A/idata/emission/smoke/GFED/GFEDV3/GFED3_3HOURLY/"
    datir          = "/Volumes/TOSHIBA_3A/idata/emission/smoke/QFED/v2.4r6/0.25/"
    datir_grid     = "../results/grids/"
    datir_out      = "../results/regrid/QFED/"
    idom           = "d01"   ;"d02"   ; d01 ; domain index
    nx             = 129     ;297     ; 129
    ny             = 84      ;120     ; 84

;   =======================================================================================
;   BECAUSE GLOBAL DATA IS VERY LARGE. SO SUBSET IS SELECTED. FOR 0.3125x0.25, SAHEL CASE
;   FROM -10N (INDEX = 320) TO 20N (INDEX = 440) FOR LATITUDE,  SO 440-320+1 = 121
;   FROM -20E (INDEX = 512) TO 45E (INDEX = 720) FOR LONGITUDE, SO 720-512+1 = 209
;   =======================================================================================
    dlon           = 0.3125
    dlat           = 0.25
    mt             = 8
    mx1            = 209
    my1            = 121
    mx2            = 720
    my2            = 360
    start_year     = 2010
    start_month    = 1
    end_year       = 2010
    end_month      = 12
    l_logic        = False
;   CONVERT kg/m^2/s TO g/m^2/hr
    factor         = 1000.*86400.
;   end modification region
    bv2            = -9999.0

    ; read in the wrf lat and lon info
    rlat = asciiread(datir_grid + "lat_"+idom+".txt", (/ny,nx/), "float")
    rlon = asciiread(datir_grid + "lon_"+idom+".txt", (/ny,nx/), "float")
    ; note: corner_rlat and corner_rlon only for "conserve" interpolation method
    corner_rlat = asciiread(datir_grid + "corner_lat_"+idom+".txt", (/ny,nx,4/), "float")
    corner_rlon = asciiread(datir_grid + "corner_lon_"+idom+".txt", (/ny,nx,4/), "float")

    aernm    = (/"oc","bc","pm25"/)
    ndate    = 365
    gfed_3hr = new((/mt,my2,mx2/),"float")

; below to begin the regridding from lat-lon to lambert projection
    srcGridName1= "src_SCRIP1.nc"
    srcGridName2= "src_SCRIP2.nc"
    dstGridName = "dst_SCRIP.nc"
    wgtFileName1= "LL_1_CWRF.nc"
    wgtFileName2= "LL_2_CWRF.nc"

    ; for source grids of QFED
    Opt                = True
    ; 0.1x0.1grid
;   Opt@LLCorner       = (/  -89.95, -179.95/)
;   Opt@URCorner       = (/   89.95,  179.95/)
    ; 0.25x0.3125grid
;   Opt@LLCorner       = (/  -90, -180.000/)
;   Opt@URCorner       = (/   90,  179.6875/)
    Opt@LLCorner       = (/  -10, -20.0/)
    Opt@URCorner       = (/   20,  45.0/)
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    latlon_to_SCRIP(srcGridName1,"0.25x0.3125",Opt)
    delete(Opt)

    ; for source grids of GFED emission fraction
    Opt                = True
    ; 0.5x0.5grid
    Opt@LLCorner       = (/  -89.75, -179.75/)
    Opt@URCorner       = (/   89.75,  179.75/)
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    latlon_to_SCRIP(srcGridName2,"0.5x0.5",Opt)
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
    Opt                = True
    Opt@InterpMethod   = "conserve"   ; "bilinear","patch","conserve"
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Debug          = True
    ESMF_regrid_gen_weights(srcGridName1, dstGridName, wgtFileName1, Opt)
    ESMF_regrid_gen_weights(srcGridName2, dstGridName, wgtFileName2, Opt)
    delete(Opt)

    ndays              = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    temp_regrid1       = new((/mt, ny, nx/), "float", bv2)
    temp_regrid2       = new((/mt, ny, nx/), "float", bv2)
    temp               = new((/my1, mx1/), "float",   bv2)

;   do iy = start_year, end_year
       cyear = iy
;      do im = start_month, end_month
          cmonth = int2string_format(im, l_logic, 0, 2)
 
          if (isleapyear(iy)) then
             ndays(1) = 29
          else
             ndays(1) = 28
          end if

          src_file1 = datir_scale3hr + cyear + "/fraction_emissions_"+cyear+cmonth+".nc"
          s1        = addfile(src_file1, "r")
          gfed_3hr  = s1->Fraction_of_Emissions

;         do id = 1, ndays(im-1)
             cday = int2string_format(id, l_logic, 0, 2)

             do iname = 0,2   ; 0,2 =>"OC", "BC", and "PM2p5"
                flnm     = datir     + "Y" + cyear + "/M" + cmonth   +  \
                           "/qfed2.emis_" + aernm(iname) + ".005." +    \
                           cyear + cmonth + cday + ".nc4"
                flnm_out = datir_out + cyear + "/" +                    \
                           "/qfed2.emis_" + aernm(iname) + ".005." +    \
                           cyear + cmonth + cday + "_" + idom + ".nc"

                print("Now Processing : " + flnm)
                ;read in the source data (this data from north to south)
                a        = addfile(flnm, "r")
                temp     = a->biomass(0, 320:440, 512:720)
                if (all(ismissing(temp))) then
                   temp = 0.0
                end if
                ; for final data
                Opt                = True
                Opt@PrintTimings   = True
                Opt@ForceOverwrite = True
                Opt@Debug          = True
                Opt@DstGridType    = "rectilinear"
                do ihr = 0, mt-1
                   temp_regrid1(ihr, :, :)= ESMF_regrid_with_weights(temp, wgtFileName1, Opt)
                   temp_regrid2(ihr, :, :)= ESMF_regrid_with_weights(gfed_3hr(ihr, :, :), wgtFileName2, Opt)
                end do

                do ihr = 0,mt-1
                   do j = 0,ny-1
                      ; unit: g/m2/3hr
                      temp_regrid1(ihr,j,:) = temp_regrid1(ihr, j, :)*temp_regrid2(ihr,j,:)*factor
                   end do
                end do

                temp_regrid1@lat2d  = rlat
                temp_regrid1@lon2d  = rlon

                system("rm -f "+flnm_out)
                n = addfile(flnm_out, "c")
                fAtt               = True            ; assign file attributes
                fAtt@CREATION_DATE = systemfunc ("date")
                fAtt@DECODE_BY     = "ESMF_regridding.ncl was used."
                fileattdef( n, fAtt )            ; copy file attributes
                ; define dimension
                dimNames = (/"time", "lat",  "lon"  /)
                dimSizes = (/ -1,       ny,   nx    /)
                dimUnlim = (/True,   False,   False /)
                filedimdef(n,dimNames,dimSizes,dimUnlim)

                ; output variables
                filevardef(n, "smoke_emission", "float", (/ "time", "lat", "lon" /) )
                temp_regrid1@long_name    = " QFED smoke PM2.5 emission"
                temp_regrid1@units        = "g/m2/3hr"
                temp_regrid1@missingValue = -9999.
                temp_regrid1@_FillValue   = -9999.
                n->smoke_emission         = temp_regrid1

                filevardef(n, "rlat", "float", (/ "lat", "lon" /) )
                rlat@long_name = "latitude"
                rlat@units     = "degree"
                n->rlat        = rlat

                filevardef(n, "rlon", "float", (/ "lat", "lon" /) )
                rlon@long_name = "longitude"
                rlon@units     = "degree"
                n->rlon        = rlon
             end do
;         end do
;      end do
;   end do

end 
