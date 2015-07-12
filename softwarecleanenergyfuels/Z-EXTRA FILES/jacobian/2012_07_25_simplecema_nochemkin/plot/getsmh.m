function smh=getsmh(T)


RU=83145100; %erg/(mol*K)
smh=gets(T)/RU - geth(T)/(RU*T);








end