# atmospheric delay correction testing
# jemez_insar_swe_2-12_2-19_HH
# HH
# 3/22

#geolocating is right, talk to HP about it

library(gridExtra)
library(data.table)
library(gdalUtils)
library(sp)
library(caTools)
library(rgdal)
library(rgeos)
library(ggplot2)
library(raster)
library(zoo)

# import i_angle raster
i_angle_deg_raw <-raster("/Volumes/JT/projects/uavsar/jemez/inc/i_angle_deg.tif")
plot(i_angle_deg_raw)
i_angle_deg_raw

#bring in UAVSAR rasters
files <-list.files("/Volumes/JT/projects/uavsar/jemez/atm_correct/", pattern = ".tiff", full.names = TRUE)
files
stack_raw <-stack(files)
stack_raw # inspect

####### set no data values to NA for all 5 layers, restack

# amp1
amp1 <-stack_raw[[1]]
values(amp1)[values(amp1) == 0] = NA
plot(amp1)
hist(amp1)

# amp2
amp2 <-stack_raw[[2]]
values(amp2)[values(amp2) == 0] = NA
plot(amp2)
hist(amp2)

# cor
cor <-stack_raw[[3]]
values(cor)[values(cor) == 0] = NA
plot(cor)
hist(cor)

# dem
dem <-stack_raw[[4]]
values(dem)[values(dem) == -10000] = NA
plot(dem)
hist(dem)

# unw
unw <-stack_raw[[5]]
values(unw)[values(unw) == 0] = NA
plot(unw)
hist(unw)

insar_stack <-stack(amp1, amp2, cor, dem, unw)
insar_stack
plot(insar_stack)



####################################

# create snow masks
{
# bring in single image that is what we need
fsca_raw <-raster("/Volumes/JT/projects/uavsar/jemez/fsca/02_18_2020/LC08_CU_010012_20200218_20200227_C01_V01_SNOW/LC08_CU_010012_20200218_20200227_C01_V01_SNOW.tif")
fsca <-projectRaster(fsca_raw,
                     crs=crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fsca #check 
fsca_raw # compare
values(fsca)[values(fsca) == 0] = NA

#test plot with DEM, this looks good but resolution all off?
plot(dem)
plot(fsca, add=TRUE)

#crop to DEM
fsca_crop <-crop(fsca, insar_stack)
plot(fsca_crop)
hist(fsca_crop)

plot(dem)
plot(fsca_crop, add=TRUE)

# resample to get down to DEM res and save it
fsca_crop_resamp <- resample(fsca_crop, insar_stack, method='bilinear')
fsca_crop_resamp

# mask for cor values and save it
fsca_final <-mask(fsca_crop_resamp, cor, maskvalue = NA)
plot(cor)
plot(fsca_final, add=TRUE)
writeRaster(fsca_final, "/Volumes/JT/projects/uavsar/jemez/atm_correct/fsca.tif")

#create snow mask for pixels that are over 90 percent snow
snow_mask <-fsca_final
values(snow_mask)[values(snow_mask) > 0] = 1
plot(snow_mask)
#writeRaster(snow_mask, "/Volumes/JT/projects/uavsar/jemez/atm_correct/snow_mask.tif")

#create no snow mask for pixels that are over 90 percent snow
no_snow_mask <-fsca_final
values(no_snow_mask)[values(is.na(no_snow_mask))] = -999
values(no_snow_mask)[values(no_snow_mask) > 0] = NA
values(no_snow_mask)[values(no_snow_mask) < 0] = 1
no_snow_mask <-mask(no_snow_mask, cor, maskvalue = NA)
plot(no_snow_mask)
hist(no_snow_mask)
plot(snow_mask, add = TRUE)
#writeRaster(no_snow_mask, "/Volumes/JT/projects/uavsar/jemez/atm_correct/no_snow_mask.tif")

}

# snow mask
no_snow_mask <-raster( "/Volumes/JT/projects/uavsar/jemez/atm_correct/no_snow_mask.tif")
snow_mask <-raster( "/Volumes/JT/projects/uavsar/jemez/atm_correct/snow_mask.tif")

#mask phase for snow and no snow
unw_snow <- mask(unw, snow_mask, maskvalue = NA)
unw_no_snow <-mask(unw, no_snow_mask, maskvalue = NA)

plot(unw_no_snow)
hist(unw_no_snow)
plot(unw_snow)
hist(unw_snow)


############################
### plot no snow phase and no phase, see if signal is from ATM
###########################

lm_snow <-lm(unwrapped_phase ~ x, unw_snow_points)
lm_no_snow <-lm(unwrapped_phase ~ x, unw_no_snow_points)

######### snow
# convert to points for graph
unw_snow_points <-as.data.frame(rasterToPoints(unw_snow))
head(unw_snow_points)
colnames(unw_snow_points)[3] <- "unwrapped_phase"
unw_snow_points <-filter(unw_snow_points, unwrapped_phase < 13)
max(unw_snow_points$unwrapped_phase)

# plot corrected data
theme_set(theme_light(base_size =12))
p9 <-ggplot(unw_snow_points, aes(x, unwrapped_phase)) +
  geom_hex(bins = 25) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  #stat_smooth_func2(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_smooth(method = "lm", se = FALSE) +
  #geom_label(
   # label="y = -20.04x - 2131.8", 
    #x=4.1,
    #y=20,
    #label.size = 0.35,
    #color = "black"
  #)+
  #geom_abline(slope = coef(lm_snow)[[2]], intercept = coef(lm_snow)[[1]], size = 1)+
  #scale_y_continuous(breaks = seq(-5,6,2))+
  labs(title = "Jemez River SCA Unwrapped Phase",
       x = "Longitude Change (deg)",
       y = "Unwrapped Phase (radians)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) + 
  annotate("text", x=-106.4, y=11, label = "y = -20.04x - 2131.8")

print(p9)

setwd("/Volumes/JT/projects/uavsar/jemez/atm_correct/")
ggsave(p9,
  file = "unw_snow_ramp.png",
  width = 6, 
  height = 4,
  dpi = 400)

######### no snow
# convert to points for graph
unw_no_snow_points <-as.data.frame(rasterToPoints(unw_no_snow))
head(unw_no_snow_points)
colnames(unw_no_snow_points)[3] <- "unwrapped_phase"
max(unw_no_snow_points$unwrapped_phase)

# plot corrected data
theme_set(theme_light(base_size =12))
p10 <-ggplot(unw_no_snow_points, aes(x, unwrapped_phase)) +
  geom_hex(bins = 25) +
  scale_fill_gradient(low = "white", high = "goldenrod") +
  #stat_smooth_func2(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_smooth(method = "lm", se = FALSE) +
  #geom_abline(slope = coef(lm_fit)[[2]], intercept = coef(lm_fit)[[1]], size = 1)+
  #scale_y_continuous(breaks = seq(-5,6,2))+
  labs(title = "Jemez River Snow Free Unwrapped Phase",
       x = "Longitude Change (deg)",
       y = "Unwrapped Phase (radians)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  annotate("text", x=-106.4, y=11, label = "y = -18.12x - 1928.4")

print(p10) #+ annotate("text", x=-106.4, y=17, label = "y = -18.12x - 1928.4")

lm_no_snow$coefficients

ggsave(p10,
       file = "unw_no_snow_ramp.png",
       width = 6, 
       height = 4,
       dpi = 400)


######################################################
#################### transect plots ##################
######################################################

# read in east west profile for phase
EW <-read.csv("/Volumes/JT/projects/uavsar/jemez/atm_correct/east-west_transect.csv")

#change 0 phase values to NA
#unw_EW$unw[unw_EW$unw == 0] <-NA
unw_10p_mean <-as.data.frame(rollmean(EW$unw, k = 10)) #calc 10 pixel rolling mean
colnames(unw_10p_mean)[1] <- "rolling_mean10" #change name
add_na <-as.data.frame(rep(NA, 9)) # add 9 NAs
colnames(add_na)[1] <- "rolling_mean10"
unw_10p_mean <-rbind(add_na, unw_10p_mean)

# 20p
unw_20p_mean <-as.data.frame(rollmean(EW$unw, k = 20)) #calc 10 pixel rolling mean
colnames(unw_20p_mean)[1] <- "rolling_mean20" #change name
add_na2 <-as.data.frame(rep(NA, 19)) # add 9 NAs
colnames(add_na2)[1] <- "rolling_mean20"
unw_20p_mean <-rbind(add_na2, unw_20p_mean)


EW_new <-cbind(EW, unw_10p_mean$rolling_mean10, unw_20p_mean$rolling_mean20)
colnames(EW_new)[6] <- "rolling_mean10"
colnames(EW_new)[7] <- "rolling_mean20"

#unw_profile <-new_df %>% arrange(desc(lon_change))
write.csv(EW_new, "/Volumes/JT/projects/uavsar/jemez/atm_correct/EW_means.csv")

##### plot ew phase data

# use 4 bc phase max and 3500 bc ele max to rescale data

yaht <-ggplot(EW_new, aes(x =lon)) +
  geom_line(aes(y=unw), size = .1) +
  geom_line(aes(y=ele * 4 / 3500), color = "darkgreen")+
  geom_line(aes(y=rolling_mean10), color = "red") + 
  geom_line(aes(y=rolling_mean20), color = "blue") +
  labs(title = "Jemez (E-W) Phase and Elevation Profile",
       x = "Longitude (deg)",
       y = "Unwrapped Phase (radians)")+
  scale_y_continuous(
    name = "Longitude Change (deg)", 
    sec.axis = sec_axis(~ .* (3500 / 4), name = "Elevation (m)"), limits = c(-2, 4))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
plot(yaht)

setwd("/Volumes/JT/projects/uavsar/jemez/atm_correct/plots")
ggsave(yaht,
  file = "ew_ele_phase_dualaxis.png",
  width = 6, 
  height = 4,
  dpi = 400)

### ew unw

ew_unw_plot <-ggplot(EW_new, aes(x =lon)) +
  geom_line(aes(y=unw), size = .1) +
  geom_line(aes(y=rolling_mean10), color = "red") + 
  geom_line(aes(y=rolling_mean20), color = "blue") +
  labs(title = "Jemez (E-W) Phase and Elevation Profile",
       x = "Longitude (deg)",
       y = "Unwrapped Phase (radians)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
plot(ew_unw_plot)


# plot ew ele data

ew_ele_plot <-ggplot(EW_new, aes(x =lon)) +
  geom_line(aes(y=ele), color = "darkgreen") +
  labs(#title = "Jemez (E-W) Elevation Profile",
       x = "Longitude (deg)",
       y = "Elevation (meters)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

ew_dual_plot <-grid.arrange(ew_unw_plot, ew_ele_plot, ncol=2)
ggsave(ew_dual_plot,
       file = "ew_ele_phase_dualplot.png",
       width = 10, 
       height = 4,
       dpi = 400)

ew_dual_plot2 <-grid.arrange(ew_unw_plot, ew_ele_plot, ncol=1)
ggsave(ew_dual_plot2,
       file = "ns_ele_phase_dualplot2.png",
       width = 10, 
       height = 5,
       dpi = 400)


#################
## north south
##################

NS <-read.csv("/Volumes/JT/projects/uavsar/jemez/atm_correct/north-south_transect.csv")

#change 0 phase values to NA
#unw_EW$unw[unw_EW$unw == 0] <-NA
unw_10p_mean <-as.data.frame(rollmean(NS$unw, k = 10)) #calc 10 pixel rolling mean
colnames(unw_10p_mean)[1] <- "rolling_mean10" #change name
add_na <-as.data.frame(rep(NA, 9)) # add 9 NAs
colnames(add_na)[1] <- "rolling_mean10"
unw_10p_mean <-rbind(add_na, unw_10p_mean)

# 20p
unw_20p_mean <-as.data.frame(rollmean(NS$unw, k = 20)) #calc 10 pixel rolling mean
colnames(unw_20p_mean)[1] <- "rolling_mean20" #change name
add_na2 <-as.data.frame(rep(NA, 19)) # add 9 NAs
colnames(add_na2)[1] <- "rolling_mean20"
unw_20p_mean <-rbind(add_na2, unw_20p_mean)


NS_new <-cbind(NS, unw_10p_mean$rolling_mean10, unw_20p_mean$rolling_mean20)
colnames(NS_new)[6] <- "rolling_mean10"
colnames(NS_new)[7] <- "rolling_mean20"




##### plot NS phase data

# use 4 bc phase max and 3500 bc ele max to rescale data

hmmt <-ggplot(NS_new, aes(x =lat)) +
  geom_line(aes(y=unw), size = .1) +
  geom_line(aes(y=ele * 4 / 3500), color = "darkgreen")+
  geom_line(aes(y=rolling_mean10), color = "red") + 
  geom_line(aes(y=rolling_mean20), color = "blue") +
  labs(title = "Jemez (N-S) Phase and Elevation Profile",
       x = "Latitude (deg)",
       y = "Unwrapped Phase (radians)")+
  scale_y_continuous(
    name = "Longitude Change (deg)", 
    sec.axis = sec_axis(~ .* (3500 / 4), name = "Elevation (m)"), limits = c(-2, 4))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
plot(hmmt)

hmmt2 <-hmmt+geom_line(aes(x =lon, y=(ele / 10)), color = "darkgreen")
plot(hmmt)

setwd("/Volumes/JT/projects/uavsar/jemez/atm_correct/")
ggsave(hmmt,
       file = "ns_ele_phase_dualaxis.png",
       width = 6, 
       height = 4,
       dpi = 400)


# plot ew ele data

ns_ele_plot <-ggplot(NS_new, aes(x =lat)) +
  geom_line(aes(y=ele), color = "darkgreen") +
  labs(#title = "Jemez (N-S) Phase and Elevation Profile",
       x = "Latitude (deg)",
       y = "Elevation (meters)")+
  ylim(c(2000,3500))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
plot(ns_ele_plot)

## unw
ns_unw_plot <-ggplot(NS_new, aes(x =lat)) +
  geom_line(aes(y=unw), size = .1) +
  geom_line(aes(y=rolling_mean10), color = "red") + 
  geom_line(aes(y=rolling_mean20), color = "blue") +
  labs(title = "Jemez (N-S) Phase and Elevation Profile",
       x = "Latitude (deg)",
       y = "Unwrapped Phase (radians)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
plot(ns_unw_plot)

ns_dual_plot <-grid.arrange(ns_unw_plot, ns_ele_plot, ncol=2)
ggsave(ns_dual_plot,
       file = "ns_ele_phase_dualplot.png",
       width = 10, 
       height = 4,
       dpi = 400)

ns_dual_plot2 <-grid.arrange(ns_unw_plot, ns_ele_plot, ncol=1)
ggsave(ns_dual_plot2,
       file = "ns_ele_phase_dualplot2.png",
       width = 10, 
       height = 5,
       dpi = 400)

#




















# extract lat lon information for multiplication
lon <- unw_zero_na
lat <- unw_zero_na
xy <- coordinates(unw_zero_na)
lon[] <- xy[, 1]
lat[] <- xy[, 2]
plot(lon)
plot(lat)
#writeRaster(lat,"/Volumes/JT/projects/uavsar/jemez/swe_calc/unw_lat.tif")
#writeRaster(lon,"/Volumes/JT/projects/uavsar/jemez/swe_calc/unw_lon.tif")

# correct for slope using best fit equation

#### lm results
# (Intercept)           x 
# -2115.1265    -19.8767 

slope_correct <-function(unw, lon){return((unw - ((lon * -19.8767) - 2115.1265)))}
unw_corrected <-slope_correct(unw_zero_na, lon)
plot(unw_corrected)
#writeRaster(unw_corrected,"/Volumes/JT/projects/uavsar/jemez/swe_calc/unw_corrected.tif")

# test hist
hist(unw_corrected)
hist(unw_zero_na)

# convert to points for graph
unw_corrected_points <-as.data.frame(rasterToPoints(unw_corrected))
head(unw_corrected_points)
colnames(unw_corrected_points)[3] <- "unwrapped_phase"

# plot corrected data
theme_set(theme_light(base_size =12))
p9 <-ggplot(unw_corrected_points, aes(x, unwrapped_phase)) +
  geom_hex(bins = 25) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  #stat_smooth_func2(geom="text",method="lm",hjust=0,parse=TRUE) +
  #geom_smooth(method = "lm", se = FALSE) +
  #geom_abline(slope = coef(lm_fit)[[2]], intercept = coef(lm_fit)[[1]], size = 1)+
  #scale_y_continuous(breaks = seq(-5,6,2))+
  labs(title = "Jemez River Unwrapped Phase vs. Longitude Corrected",
       x = "Longitude Change (deg)",
       y = "Unwrapped Phase (radians)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
print(p9)

setwd("/Volumes/JT/projects/uavsar/jemez/swe_calc/")
#ggsave(p9,
file = "unw_corrected_graph.png",
width = 6, 
height = 4,
dpi = 400)



########################################################
### calculating 
########################################################

density <- .29 # get a real number and do senativity analysis
di_elc <- 1.4 # 
wL <- 23.8403545

# first step
insar_constant <-function(inc){((-4*pi)/wL)*(cos(inc) - sqrt(di_elc - sin((inc)^2)))}
insar_constant_rast <-insar_constant(i_angle_deg)
hist(insar_constant_rast)
plot(insar_constant_rast)

#do swe change calc
delta_swe_rast <-insar_constant_rast*unw_corrected
plot(delta_swe_rast)
hist(delta_swe_rast, breaks = 100)
#writeRaster(delta_swe_rast,"/Volumes/JT/projects/uavsar/jemez/swe_calc/delta_swe_raster.tif")

#mask for snow
delta_swe_snow_mask <- mask(delta_swe_rast, snow_mask, maskvalue = NA)
plot(delta_swe_snow_mask)
#writeRaster(delta_swe_snow_mask,"/Volumes/JT/projects/uavsar/jemez/swe_calc/delta_swe_snow_mask.tif")
hist(delta_swe_snow_mask)
freq(delta_swe_snow_mask, digits =0)

cc <-raster("/Volumes/JT/projects/uavsar/jemez/nlcd/cc_final.tif")
hist(cc)
values(cc)[values(cc) > 1] = NA
cc_swe_mask <- mask(delta_swe_snow_mask, cc, maskvalue = NA)
plot(cc_swe_mask)
#writeRaster(cc_swe_mask, "/Volumes/JT/projects/uavsar/jemez/swe_calc/cc_swe_mask.tif")


####################################### 
###### find no change point and subtract
#######################################

#import
delta_swe_rast <-raster("/Volumes/JT/projects/uavsar/jemez/swe_calc/delta_swe_raster.tif")

#define no change point phase
no_change_point_phase <-0.268596

#subtract from raster
delta_swe_abs <-delta_swe_rast-no_change_point_phase
plot(delta_swe_abs)

#mask for snow cover
delta_swe_snow_mask_abs <- mask(delta_swe_abs, snow_mask, maskvalue = NA)
plot(delta_swe_snow_mask_abs)
writeRaster(delta_swe_snow_mask_abs,"/Volumes/JT/projects/uavsar/jemez/swe_calc/delta_swe_snow_mask_abs.tif")

delta_swe_snow_mask_abs1 <-raster("/Volumes/JT/projects/uavsar/jemez/swe_calc/2-12_2-19/delta_swe_snow_mask_abs_2-12_2-19_HH.tif")

hist(delta_swe_snow_mask_abs1,
     breaks = 100,
     main= "Jemez Change in SWE 2/12-2/19 HH",
     xlab= "SWE Change (cm)",
     xlim=c(-3,2),
     col= "darkgreen"
)
