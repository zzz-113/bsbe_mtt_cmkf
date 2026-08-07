function soundSpeed_t = getWBAT_SoundSpeed(t)

dat1 = load('WBAT_CTD_2021.mat');
dat2 = load('WBAT_CTD_2022.mat');

time = datetime([dat1.Dat.mtime dat2.Dat.mtime],'ConvertFrom','datenum');
temp = [dat1.Dat.temp dat2.Dat.temp];
salt = [dat1.Dat.salt dat2.Dat.salt];
pres = [dat1.Dat.press dat2.Dat.press];

cw = gsw_sound_speed(salt,temp,pres);


soundSpeed_t = interp1(time,cw,t);


end