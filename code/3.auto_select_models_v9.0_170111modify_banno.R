
# adjust memory -------------------------------------------------

rm(list=ls())
sapply(X=1:9, FUN=gc)
# memory.size(max=T) ## codes in Windows
# history(Inf)

# setting -----------------------------------------------------------------

SET_RESPOSITORY_PATH=
        "/Users/jie.i.zhang/_eny/dili/&&proj/16.07-08 Kubota/code_examine_170110"
subFd="auto_model_select_v9.0"

SET_CURRENT_MONTH="2016/06"

frcst_mon_future=12
frcst_mon_testing=3
# frcst_mon_testing=min(12,nrow(xy.sub)-12)

# path ----------------------------------------------------------

SET_PATH=list(
        SRC_DATA=paste0(SET_RESPOSITORY_PATH, "/src_9.0/"),
        ASSIST=paste0(SET_RESPOSITORY_PATH, "/assist/"),
        ALL_RLT=paste0(SET_RESPOSITORY_PATH, "/rlt/"),
        RLT=paste0(SET_RESPOSITORY_PATH, "/rlt/",subFd),
        MID_RLT=paste0(SET_RESPOSITORY_PATH, "/rlt/",subFd,"/mid/"),
        PLOT=paste0(SET_RESPOSITORY_PATH, "/rlt/",subFd,"/plot/")
)
setwd(SET_RESPOSITORY_PATH)

if (!dir.exists(SET_PATH$ALL_RLT)) {dir.create(SET_PATH$ALL_RLT)}
if (!dir.exists(SET_PATH$RLT)) {dir.create(SET_PATH$RLT)}
if (!dir.exists(SET_PATH$MID_RLT)) {dir.create(SET_PATH$MID_RLT)}
if (!dir.exists(SET_PATH$PLOT)) {dir.create(SET_PATH$PLOT)}

# library -------------------------------------------------------

SET_PACKAGE_LIST=c(
        "dplyr","tidyr",
        "car","zoo",
        "readr",
        "psych","forecast",
        "e1071","gbm","randomForest","RSNNS",
        # "foreach",
        "ggplot2"
)
to_install_packages=setdiff(SET_PACKAGE_LIST,installed.packages()[,1])
if (length(to_install_packages)!=0) {
        install.packages(to_install_packages)
}
if (all(sapply(
        X=SET_PACKAGE_LIST,
        FUN=require,
        character.only=T,
        warn.conflicts=F,
        quietly=T
))){print('load packages correct')} else{
        print('load packages wrong')     
}

# func ----------------------------------------------------------

mon_minus <- function(t_end,t_start) {
        t_end_yr=as.integer(strsplit(t_end,"/")[[1]][1])
        t_end_mn=as.integer(strsplit(t_end,"/")[[1]][2])
        t_start_yr=as.integer(strsplit(t_start,"/")[[1]][1])
        t_start_mn=as.integer(strsplit(t_start,"/")[[1]][2])
        mon_m=t_end_mn-t_start_mn + (t_end_yr-t_start_yr)*12 + 1
        return(mon_m)
}

write__csv <- function(outD,outN){
        outP=paste0(SET_PATH$MID_RLT,outN,".csv")
        write.csv(
                outD,
                outP,
                row.names=F,
                na="NA",
                fileEncoding = "SHIFT-JIS"
        )
}

# modelling -----------------------------------------------------

###############################################################
## arima

am_func=function(tr_data,f_mon){
        tr_data.ts=ts(tr_data,frequency = 12)
        frcst_value=try(
                forecast.Arima(
                        auto.arima(
                                tr_data.ts,
                                lambda=NULL #lambda
                        ),
                        h=f_mon,
                        lambda=NULL
                )$mean
        )
        if(class(frcst_value)=="try-error") {
                frcst_value=rep(NA,f_mon)
        }
        return(as.numeric(frcst_value))
}

###############################################################
## random forest

rf_func=function(tr_data,f_mon){
        tr_data=as.data.frame(tr_data)
        set.seed(2189)
        mtry_fc <- try(
                tuneRF(
                        x=tr_data[,
                                  which(names(tr_data)=="x1"):
                                          which(names(tr_data)=="x12")
                                  ],
                        y=tr_data$qnty_normal,
                        # mtryStart=1,
                        ntreeTry=200,
                        stepFactor=1.5,
                        improve=0.01,
                        trace=F,
                        plot=F
                )
                ,silent = T)
        if (class(mtry_fc)!="try-error") {
                frcst_value=rep(NA,f_mon)
                best.m_fc <- mtry_fc[which.min(mtry_fc[, 2]), 1]
                for (i in 1:f_mon) {
                        last_row=nrow(tr_data)
                        tr_data[last_row+1,]=c(
                                tr_data$cd[1],
                                tr_data$clst_type[1],
                                tr_data$yr_mon[1],
                                tr_data$qnty_min[1],
                                tr_data$qnty_max[1],
                                NA,
                                tr_data[last_row,
                                        which(names(tr_data)=="qnty_normal"):
                                                which(names(tr_data)=="x11")
                                        ]
                        )
                        set.seed(2189)
                        frcst_value[i]=predict(
                                randomForest(
                                        qnty_normal~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12,
                                        data=tr_data[-nrow(tr_data),],
                                        mtry=best.m_fc
                                )
                                ,newdata=tr_data[(last_row+1),
                                                 which(names(tr_data)=="x1"):
                                                         which(names(tr_data)=="x12")
                                                 ]
                        )
                        tr_data$qnty_normal[nrow(tr_data)]=frcst_value[i]
                }
        } else {
                return(rep(NA,f_mon))
        }
        return(frcst_value)
}

###############################################################
## gradient boosting


        frcst_test=frcst_test %>% mutate(
                gbm=gbm_func(
                        data.frame(xy.sub[c(1:(nrow(xy.sub)-frcst_mon_testing)),])
                        ,frcst_mon_testing
                )
        )

gbm_func=function(tr_data,f_mon){
        set.seed(2189)
        fit.gbm_fc=try(
                gbm(
                        qnty_normal~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12,
                        # qnty_normal~x1,
                        data=tr_data[,
                                     which(names(tr_data)=="qnty_normal"):
                                             which(names(tr_data)=="x12")
                                     ],
                        distribution="gaussian",
                        n.trees=5000, # the number of trees B
                        shrinkage=0.001, # shrinkage parameter lambda
                        interaction.depth=2, # the number d of splits in each tree
                        cv.folds=5,
                        bag.fraction = 0.7,
                        n.minobsinnode=1
                        
                )
                ,silent = T)
        # gbm.perf(fit.gbm_fc)
        if (class(fit.gbm_fc)!="try-error") {
                frcst_value=rep(NA,f_mon)
                for (i in 1:f_mon) {
                        last_row=nrow(tr_data)
                        tr_data[last_row+1,]=c(
                                tr_data$cd[1],
                                tr_data$clst_type[1],
                                tr_data$yr_mon[1],
                                tr_data$qnty_min[1],
                                tr_data$qnty_max[1],
                                NA,
                                tr_data[last_row,
                                        which(names(tr_data)=="qnty_normal"):
                                                which(names(tr_data)=="x11")
                                        ]
                        )
                        set.seed(2189)
                        frcst_value[i]=predict(
                                fit.gbm_fc
                                ,newdata=tr_data[(last_row+1),
                                                 which(names(tr_data)=="x1"):
                                                         which(names(tr_data)=="x12")
                                                 ]
                                ,n.trees=gbm.perf(fit.gbm_fc,plot.it=F)
                        )
                        tr_data$qnty_normal[nrow(tr_data)]=frcst_value[i]
                }
        } else {
                return(rep(NA,f_mon))
        }
        return(frcst_value)
}

###############################################################
## Support Vector Regression

svr_func=function(tr_data,f_mon){
        set.seed(2189)
        tuneResult_fc= try(
                tune(
                        svm,
                        qnty_normal~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12,
                        data=tr_data[,
                                     which(names(tr_data)=="qnty_normal"):
                                             which(names(tr_data)=="x12")
                                     ],
                        ranges = list(
                                epsilon = seq(0,1,0.01)
                                # ,cost = 2^(2:9)
                        )
                )
                ,silent = T)
        if (class(tuneResult_fc)!="try-error") {
                frcst_value=rep(NA,f_mon)
                for (i in 1:f_mon) {
                        last_row=nrow(tr_data)
                        tr_data[last_row+1,]=c(
                                tr_data$cd[1],
                                tr_data$clst_type[1],
                                tr_data$yr_mon[1],
                                tr_data$qnty_min[1],
                                tr_data$qnty_max[1],
                                NA,
                                tr_data[last_row,
                                        which(names(tr_data)=="qnty_normal"):
                                                which(names(tr_data)=="x11")
                                        ]
                        )
                        set.seed(2189)
                        frcst_value[i]=predict(
                                tuneResult_fc$best.model
                                ,newdata=tr_data[(last_row+1),
                                                 which(names(tr_data)=="x1"):
                                                         which(names(tr_data)=="x12")
                                                 ]
                        )
                        tr_data$qnty_normal[nrow(tr_data)]=frcst_value[i]
                }
        } else {
                return(rep(NA,f_mon))
        }
        return(frcst_value)
}

###############################################################
## Elman NN

enn_func=function(tr_data,f_mon){
        set.seed(2189)
        fit.elman=try(
                elman(
                        tr_data[,
                                which(names(tr_data)=="x1"):
                                        which(names(tr_data)=="x12")
                                ],
                        tr_data[,
                                which(names(tr_data)=="qnty_normal")
                                ],
                        size=c(1,1),
                        learnFuncParams = c(0.01),
                        maxit = 30000
                )
        )
        # plotIterativeError(fit.elman)
        if (!is.element("try-error",class(fit.elman))) {
                frcst_value=rep(NA,f_mon)
                for (i in 1:f_mon) {
                        last_row=nrow(tr_data)
                        tr_data[last_row+1,]=c(
                                tr_data$cd[1],
                                tr_data$clst_type[1],
                                tr_data$yr_mon[1],
                                tr_data$qnty_min[1],
                                tr_data$qnty_max[1],
                                NA,
                                tr_data[last_row,
                                        which(names(tr_data)=="qnty_normal"):
                                                which(names(tr_data)=="x11")
                                        ]
                        )
                        set.seed(2189)
                        frcst_value[i]=predict(
                                fit.elman
                                ,newdata=tr_data[(last_row+1),
                                                 which(names(tr_data)=="x1"):
                                                         which(names(tr_data)=="x12")
                                                 ]
                        )
                        tr_data$qnty_normal[nrow(tr_data)]=frcst_value[i]
                }
        } else {
                return(rep(NA,f_mon))
        }
        return(frcst_value)
}

###############################################################
## Jordan NN

jnn_func=function(tr_data,f_mon){
        set.seed(2189)
        fit.jordan=try(
                jordan(
                        tr_data[,
                                which(names(tr_data)=="x1"):
                                        which(names(tr_data)=="x12")
                                ],
                        tr_data[,
                                which(names(tr_data)=="qnty_normal")
                                ],
                        size=2,
                        learnFuncParams = c(0.01),
                        maxit = 30000
                )
        )
        # plotIterativeError(fit.jordan)
        if (!is.element("try-error",class(fit.jordan))) {
                frcst_value=rep(NA,f_mon)
                for (i in 1:f_mon) {
                        last_row=nrow(tr_data)
                        tr_data[last_row+1,]=c(
                                tr_data$cd[1],
                                tr_data$clst_type[1],
                                tr_data$yr_mon[1],
                                tr_data$qnty_min[1],
                                tr_data$qnty_max[1],
                                NA,
                                tr_data[last_row,
                                        which(names(tr_data)=="qnty_normal"):
                                                which(names(tr_data)=="x11")
                                        ]
                        )
                        set.seed(2189)
                        frcst_value[i]=predict(
                                fit.jordan
                                ,newdata=tr_data[(last_row+1),
                                                 which(names(tr_data)=="x1"):
                                                         which(names(tr_data)=="x12")
                                                 ]
                        )
                        tr_data$qnty_normal[nrow(tr_data)]=frcst_value[i]
                }
        } else {
                return(rep(NA,f_mon))
        }
        return(frcst_value)
}


# src data ------------------------------------------------------

src_list <- as.matrix(lapply(
        X=paste0(
                SET_PATH$SRC_DATA,
                list.files(SET_PATH$SRC_DATA)
        ),
        function(x) 
                read.csv(
                        file = x,
                        header = T,
                        na.strings = "",
                        colClasses = "character",
                        encoding = "SHIFT-JIS"
                        )
))
src_data=data.frame()
for (i in 1:length(src_list)) {
        src_data=rbind(src_data,src_list[[i]])
}
## change name to English to avoid errors
colnames(src_data)=c("cd","nam","comp","yr_mon","qnty")
src_data$qnty=as.integer(src_data$qnty)
src_data=src_data %>% filter(yr_mon<=SET_CURRENT_MONTH)

## mid output (could be commented)
src_data=src_data %>% arrange(cd, yr_mon)
write__csv(src_data,"01_src_data")

# summarize -----------------------------------------------------

src.summary.all <- group_by(src_data,cd) %>% summarise(
        sal_start=first(yr_mon)
        ,sal_end=last(yr_mon)
        ,total_month=mon_minus(sal_end,sal_start)
        ,data_month=n()
        ,tmp_2016=mon_minus(SET_CURRENT_MONTH,sal_end)
        ,qnty.all=sum(qnty)
        ,qnty.mon_ave=mean(qnty,na.rm=T)
        ,qnty.mon_median=median(qnty,na.rm=T)
        ,qnty.mon_max=max(qnty,na.rm=T)
        ,qnty.mon_min=min(qnty,na.rm=T)
) %>% mutate(
        data_2016=ifelse(tmp_2016<=1,TRUE,FALSE)
        ,NA_month=total_month-data_month
        ,NA_rate=NA_month/total_month
) %>% select(
        cd,
        sal_start,
        sal_end,
        total_month,
        data_month,
        NA_month,
        NA_rate,
        data_2016,
        qnty.mon_ave,
        qnty.mon_median,
        qnty.mon_min,
        qnty.mon_max,
        qnty.all
) %>% arrange(
        cd
)
write__csv(src.summary.all,"02_date_summary")

# rule out ------------------------------------------------------

yr16_data=src_data %>% filter(
        yr_mon>="2016/01"
)
yr16_data=group_by(yr16_data,cd) %>% summarise(
        qnty_mean_yr16=mean(qnty,na.rm=T)
)

src.summary.all=src.summary.all %>% left_join(
        yr16_data,by="cd"
)
src.summary=src.summary.all %>% filter(
        data_month>=25
        ,sal_end>=SET_CURRENT_MONTH
        ,qnty.mon_ave>20,
        # ,qnty.mon_median>20
        qnty_mean_yr16>=20
)

sel_data=src.summary %>% inner_join(
        src_data %>%
                select(cd,yr_mon,qnty)
        ,by="cd"
)
sel_data=sel_data %>% left_join(
        group_by(sel_data,cd) %>% summarize(
                qnty_max=max(qnty),
                qnty_min=min(qnty)
        )
        ,by="cd"
) %>% mutate(
        qnty_normal=(qnty-qnty_min)/(qnty_max-qnty_min)
)

write__csv(sel_data,"03_sel_data")

# hierarchical cluster ------------------------------------------

# scale
xs=sel_data %>% select(
        cd,yr_mon,qnty
) %>% spread(
        yr_mon,qnty
)
xs.scale=as.data.frame(t(scale(t(xs[-1]),center = T,scale = T)))
xs.scale$cd=xs$cd

## calculate distance
hc.fit=hclust(
        dist(
                xs.scale[,-which(names(xs.scale)=="cd")]
                ,method="euclidean"
        )
        ,method="ward.D"
)

## cut cluster
nclst=5

# # plot for cluster
# plot(hc.fit)
# rect.hclust(hc.fit,k=nclst,border="red")

clst_grps=data.frame(
        seq=c(nrow(xs):1),
        cd=xs$cd,
        clst_num=cutree(hc.fit, k=nclst)
) %>% inner_join(
        read_csv(paste0(SET_PATH$ASSIST,"clst_map.csv"))
        ,by="clst_num"
)
clst_grps$cd=as.character(clst_grps$cd)
write__csv(clst_grps,"04_clst_grps")

sel_data_clst=clst_grps %>% select(
        seq,cd,clst_type
) %>% inner_join(
        sel_data %>% select(
                cd, yr_mon, qnty, qnty_max, qnty_min, qnty_normal
        )
        ,by="cd"
)

y.spread= sel_data_clst %>% select(
        cd, yr_mon, qnty_normal
) %>% spread(
        # yr_mon,qnty_normal,fill=0
        yr_mon,qnty_normal
)
write__csv(y.spread,"05_y.spread")

y.gather=y.spread %>% gather(
        yr_mon,qnty_normal,2:ncol(y.spread)
) %>% arrange(
        cd,yr_mon
) %>% left_join(
        group_by(sel_data_clst,cd) %>% summarise(
                clst_type=first(clst_type)
                ,qnty_max=first(qnty_max)
                ,qnty_min=first(qnty_min)
        ),by="cd"
)
write__csv(y.gather,"04_y.gather")

# test, select & forecast ---------------------------------------

sel_cd=sort(unique(y.spread$cd),decreasing=T)

# sel_cd=c(
# 	"T1880A2-42",
# 	"L3301HST",
# 	"80_BX25DLB-T-1"
# "L4701HST"
# ,"ZG227LA-60"
# 	,"ZG327PA-60"
# ,"ZG222A-48"
# ,"Z725KH-60"
# ,"L3560DT"
# ,"L3800HST"
# ,"B2650HSD"
# ,"L3560GST"
# ,"ZG327PA-60"
# ,"M59"
# )
# sel_cd=c(
#         "U17VR1",
#         "T1880A2-42",
#         "SVL75-2HWC",
#         "RTV1140CPX-H",
#         "RTV-X1100CWL-HS",
#         "M7060HDC12",
#         "M7060HD",
#         "L6060HSTC",
#         "L3301HST",
#         "BX2670RV60-1",
#         "BX25DLB-T-1",
#         "B2650HSD"
# )


frcst_test_all=data.frame()
frcst_test_error_all=data.frame()
frcst_future_all=data.frame()
frcst_actual_all=data.frame()
cd_collect=as.vector(NA)
cd_jump=as.vector(NA)

for (sel_cd_no in 1:length(sel_cd)) {
        
        # sel_cd_no=1
        ch_cd=sel_cd[sel_cd_no]
        print("###################################################")
        print(paste0(
                "Calculating for hb: ",
                ch_cd,
                " ( No. ",
                sel_cd_no,
                " / ",
                length(sel_cd),
                " )"
        ))
        
        y.sub=y.gather %>% filter(cd==ch_cd) %>%
                select(cd,clst_type,
                       yr_mon,qnty_min,qnty_max,qnty_normal)
        y.sub$qnty_normal=as.zoo(y.sub$qnty_normal)
        # plot(y.sub$qnty_normal,type="b")
        
        xy.sub=y.sub %>% mutate(
                x1=lag(qnty_normal,1),
                x2=lag(qnty_normal,2),
                x3=lag(qnty_normal,3),
                x4=lag(qnty_normal,4),
                x5=lag(qnty_normal,5),
                x6=lag(qnty_normal,6),
                x7=lag(qnty_normal,7),
                x8=lag(qnty_normal,8),
                x9=lag(qnty_normal,9),
                x10=lag(qnty_normal,10),
                x11=lag(qnty_normal,11),
                x12=lag(qnty_normal,12)
        )
        # nrow(xy.sub)
        
        xy.sub=na.omit(xy.sub)
        # nrow(xy.sub)
        
        if (nrow(xy.sub)<12) {
                last_num=length(cd_jump)
                cd_jump[last_num+1]=xy.sub$cd[1]
                print(paste0(
                        "____",
                        "The data amount is insufficient",
                        "-> jump calculating for ",
                        xy.sub$cd[1]
                ))
                next
        }
        
        print(paste0("____STEP 1: selecting the best model ......"))
        
        # collect valid cd
        last_num=length(cd_collect)
        cd_collect[last_num+1]=xy.sub$cd[1]
        
        train_set=c(1:(nrow(xy.sub)-frcst_mon_testing))
        test_set=c((nrow(xy.sub)-frcst_mon_testing+1):nrow(xy.sub))
        
        frcst_test=y.sub[(nrow(y.sub)-frcst_mon_testing+1):nrow(y.sub),]
        
        ###############################################################
        ## testing arima
        
        frcst_test=frcst_test %>% mutate(
                arima=am_func(
                        y.sub[c(1:(nrow(y.sub)-frcst_mon_testing)),"qnty_normal"]
                        ,frcst_mon_testing
                )
        )
        
        ###############################################################
        ## testing random forest
        
        frcst_test=frcst_test %>% mutate(
                rf=rf_func(
                        xy.sub[c(1:(nrow(xy.sub)-frcst_mon_testing)),]
                        ,frcst_mon_testing
                )
        )
        
        ###############################################################
        ## testing gradient boosting
        
        frcst_test=frcst_test %>% mutate(
                gbm=gbm_func(
                        data.frame(xy.sub[c(1:(nrow(xy.sub)-frcst_mon_testing)),])
                        ,frcst_mon_testing
                )
        )
        
        ###############################################################
        ## testing Support Vector Regression
        
        frcst_test=frcst_test %>% mutate(
                svr=svr_func(
                        data.frame(xy.sub[c(1:(nrow(xy.sub)-frcst_mon_testing)),])
                        ,frcst_mon_testing
                )
        )
        
        ###############################################################
        ##  testing Elman NN
        
        frcst_test=frcst_test %>% mutate(
                enn=enn_func(
                        data.frame(xy.sub[c(1:(nrow(xy.sub)-frcst_mon_testing)),])
                        ,frcst_mon_testing
                )
        )
        
        ###############################################################
        ## testing Jordan NN
        
        frcst_test=frcst_test %>% mutate(
                jnn=jnn_func(
                        data.frame(xy.sub[c(1:(nrow(xy.sub)-frcst_mon_testing)),])
                        ,frcst_mon_testing
                )
        )
        
        ###############################################################
        ## de-normalize
        
        for (i in
             which(names(frcst_test)=="qnty_normal"):
             which(names(frcst_test)=="jnn")
        ) {
                frcst_test[,i]=frcst_test[,i]*
                        (frcst_test$qnty_max[1]-frcst_test$qnty_min[1])+
                        frcst_test$qnty_min[1]
        }
        colnames(frcst_test)[names(frcst_test)=="qnty_normal"]="qnty"
        frcst_test_all=rbind(frcst_test_all,frcst_test)
        
        ###############################################################
        ## calculate testing metrics
        
        frcst_test_error=frcst_test %>% mutate(
                # arima_rr=abs(arima-qnty)/qnty
                # ,arima_bc_rr=abs(arima_bc-qnty)/qnty
                # ,rf_rr=abs(rf-qnty)/qnty
                # ,gbm_rr=abs(gbm-qnty)/qnty
                # ,svr_rr=abs(svr-qnty)/qnty
                # ,enn_rr=abs(enn-qnty)/qnty
                # ,jnn_rr=abs(jnn-qnty)/qnty
                arima_rss=(arima-qnty)^2
                # ,arima_bc_rss=(arima_bc-qnty)^2
                ,rf_rss=(rf-qnty)^2
                ,gbm_rss=(gbm-qnty)^2
                ,svr_rss=(svr-qnty)^2
                ,enn_rss=(enn-qnty)^2
                ,jnn_rss=(jnn-qnty)^2
        )
        frcst_test_error=group_by(frcst_test_error,cd) %>% summarise(
                clst_type=first(clst_type),
                arima_rse=sqrt(sum(arima_rss)/(n()-1)),
                # arima_bc_mse=sqrt(sum(arima_bc_rss)/(n()-1)),
                rf_rse=sqrt(sum(rf_rss)/(n())),
                gbm_rse=sqrt(sum(gbm_rss)/(n())),
                svr_rse=sqrt(sum(svr_rss)/(n())),
                enn_rse=sqrt(sum(enn_rss)/(n())),
                jnn_rse=sqrt(sum(jnn_rss)/(n()))
        ) %>% mutate(
                min_rse=min(
                        ifelse(is.na(arima_rse),Inf,arima_rse),
                        # ifelse(is.na(arima_bc_mse),999,arima_bc_mape),
                        ifelse(is.na(rf_rse),Inf,rf_rse),
                        ifelse(is.na(gbm_rse),Inf,gbm_rse),
                        ifelse(is.na(svr_rse),Inf,svr_rse),
                        ifelse(is.na(enn_rse),Inf,enn_rse),
                        ifelse(is.na(jnn_rse),Inf,jnn_rse)
                )
        )
        frcst_test_error=frcst_test_error %>% mutate(
                best_model= gsub(
                        "_rse"
                        ,""
                        ,names(frcst_test_error)[
                                which.min(
                                        frcst_test_error[1,
                                                         which(names(frcst_test_error)
                                                               =="arima_rse"):
                                                                 which(names(frcst_test_error)
                                                                       =="jnn_rse")
                                                         ]
                                )+2
                                ]
                )
        )
        print("testing result:")
        print(as.data.frame(frcst_test_error))
        
        frcst_test_error_all=rbind(frcst_test_error_all,frcst_test_error)
        
        ###############################################################
        # forecast future sales
        
        print(paste0(
                "____STEP 2: forecasting future ",
                frcst_mon_future,
                " month(s)......"
        ))
        frcst_future=as.data.frame(matrix(
                nrow=frcst_mon_future
        ))
        frcst_future[,1]=frcst_test$cd[1]
        colnames(frcst_future)="cd"
        last_year=strsplit(frcst_test$yr_mon[nrow(frcst_test)],"/")[[1]][1]
        last_month=strsplit(frcst_test$yr_mon[nrow(frcst_test)],"/")[[1]][2]
        for (j in 1:frcst_mon_future) {
                last_month=as.integer(last_month)+1
                if (last_month<10) {
                        last_month=paste0("0",last_month)
                }
                if (as.integer(last_month)>12) {
                        last_year=as.integer(last_year)+1
                        last_month="01"
                }
                frcst_future$yr_mon[j]=paste(last_year,last_month,sep="/")
        }
        frcst_future$clst_type=frcst_test$clst_type[1]
        frcst_future$qnty_min=frcst_test$qnty_min[1]
        frcst_future$qnty_max=frcst_test$qnty_max[1]
        frcst_future$best_model=frcst_test_error$best_model
        
        #  forecast future sales
        if (frcst_test_error$best_model=="-") {next}
        if (frcst_test_error$best_model=="arima"){
                am_t=am_func(
                        y.sub[,"qnty_normal"]
                        ,frcst_mon_future
                )
                for (iii in 1:nrow(frcst_future)) {
                        frcst_future$frcst_from_now[iii]=ifelse(
                                any(is.na(am_t[iii]))
                                ,
                                am_func(
                                        y.sub[c(1:(nrow(y.sub)-frcst_mon_testing)),"qnty_normal"]
                                        ,frcst_mon_testing+frcst_mon_future
                                )[-c(1:3)][iii],
                                am_t[iii]
                        )
                }
        }
        if (frcst_test_error$best_model=="rf" |
            frcst_test_error$best_model=="gbm" |
            frcst_test_error$best_model=="svr" |
            frcst_test_error$best_model=="enn" |
            frcst_test_error$best_model=="jnn"
        ) {
                eval(
                        parse(
                                text =
                                        paste0(
                                                "frcst_future$frcst_from_now=",
                                                frcst_test_error$best_model,
                                                "_func(",
                                                "data.frame(xy.sub)",
                                                ",frcst_mon_future",
                                                ") "
                                        )
                        )
                )
        }
        frcst_future$frcst_from_now=frcst_future$frcst_from_now*
                (frcst_future$qnty_max[1]-frcst_future$qnty_min[1])+
                frcst_future$qnty_min[1]
        frcst_future_all=rbind(frcst_future_all,frcst_future)
        
        frcst_future.simple=
                frcst_future[,c("cd","yr_mon",
                                "clst_type","frcst_from_now")]
        print(frcst_future.simple)
        
        frcst_actual=sel_data_clst %>% filter(
                cd==ch_cd
        ) %>% select(
                cd,yr_mon,clst_type,qnty
        ) %>% full_join(
                frcst_future.simple
                ,by=c("cd","yr_mon","clst_type")
        ) %>% arrange(
                cd,yr_mon
        )
        
        frcst_actual[
                frcst_actual$yr_mon==SET_CURRENT_MONTH
                ,]$frcst_from_now=
                frcst_actual[
                        frcst_actual$yr_mon==SET_CURRENT_MONTH
                        ,]$qnty
        
        frcst_actual_all=rbind(frcst_actual_all,frcst_actual)
        
        jpeg(
                filename = paste0(SET_PATH$PLOT,
                                  sel_cd_no,
                                  "_",
                                  y.sub$cd[1],
                                  ".jpg"),
                width=800,height=600
        )
        g <- ggplot() +
                geom_line(data = frcst_actual,
                          aes(x = as.Date(as.yearmon(
                                  as.character(frcst_actual$yr_mon), "%Y/%m"
                          ))
                          , y = qnty),
                          na.rm = T,size = 1) +
                geom_line(data = frcst_actual,
                          aes(x = as.Date(as.yearmon(
                                  as.character(frcst_actual$yr_mon), "%Y/%m"
                          ))
                          , y = frcst_from_now
                          , color = "red"),
                          na.rm = T,size = 2
                ) +
                theme(legend.position="none") +
                xlab('year/month') +
                ylab('sales quantity') +
                ggtitle(paste0(
                        "forecast ",
                        frcst_mon_future,
                        " month(s) for ",
                        frcst_actual$cd[1],
                        " (best model: ",
                        frcst_future$best_model[1],
                        ")"
                ))
        plot(g)
        dev.off()
        
}

# output --------------------------------------------------------

## output results
frcst_test_all$qnty=as.integer(frcst_test_all$qnty)

## add cluster results
if (nrow(frcst_test_all)>0) {
        write__csv(frcst_test_all,"frcst_rlt_test_period_all")
}

if (nrow(frcst_test_error_all)>0) {
        write__csv(frcst_test_error_all,"best_model_minimal_RMSE")
}

if (nrow(frcst_future_all)>0) {
        write__csv(frcst_future_all,"frcst_rlt_future_all")
}

if (nrow(frcst_actual_all)>0) {
        write__csv(frcst_actual_all,"frcst_actual_all")
}
