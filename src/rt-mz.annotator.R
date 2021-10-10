options(stringsAsFactors = FALSE)
#dir = dirname(sys.frame(1)$ofile)
dir = getwd()
atom.weights = read.table(paste0(dir,'/src/atom.weights.tab'),sep='\t',header=T)
atom.weights[atom.weights$Symbol=='H(2)','Symbol'] = 'D(2)'
atom.weights$sym = sapply(strsplit(atom.weights$Symbol,'(',TRUE),'[',1)
atom.weights = do.call(rbind,lapply(split(atom.weights,atom.weights$sym),function(x){x[order(-x$Abund)[1],]}))
LMDB = readRDS(paste0(dir,'/src/lm1_2021.Rdata'))

findMS2 = function(scans,rt=NULL,mz=NULL,drt=10,dmz=0.5,peak=NULL){
  if(!is.null(peak)){
    rt = peak$rt
    mz = peak$mz
  }
  scans[scans$retentionTime>rt-drt & scans$retentionTime<rt+drt & scans$precursorMZ<mz+dmz & scans$precursorMZ>mz-dmz,hfs]
}

findIsotopesWithCAMERA = function(x,peak.method='maxint'){
  require(CAMERA)
  ints = peakTable(x,method=peak.method)
  rownames(ints) = paste(ints$rt,ints$mz)
  
  an_xset   =  xsAnnotate(x)
  an_xset  =  groupFWHM(an_xset)
  an_xset    = findIsotopes(an_xset)
  r =  getPeaklist(an_xset)

  rownames(r) = paste(r$rt,r$mz)
  r = r[rownames(ints),]
  inx = c(1:8,(ncol(r)-2):ncol(r))
  r = list(peaks=r[,inx],intensity=as.matrix(ints[,-inx]))
  class(r) = c('sajr','list')
  r$peaks$peak.id = rownames(r$peaks) = rownames(r$intensity) = paste0('rt',round(r$peaks$rt),'mz',round(r$peaks$mz),'_',1:nrow(r$peaks))

  iso = do.call(rbind,lapply(strsplit(r$peaks$isotopes,'\\[|\\]'),function(x){
    if(length(x)==0) return(data.frame(iso.id=NA,iso.n=NA,iso.charge=NA))
    data.frame(iso.id=as.numeric(x[2]),iso.n=x[4],iso.charge=x[5])
  }))
  
  iso$dmz = iso$drt = iso$cor = NA
  iso$iso.n = c('M'=0,setNames(1:20,paste0('M+',1:20)))[iso$iso.n]
  iso$iso.charge = c('+'=1,setNames(2:10,paste0(2:10,'+')))[iso$iso.charge]
  iso$id = rownames(r$peaks)
  
  iso = split(iso,iso$iso.id)
  iso = lapply(iso,function(x){
    x = x[order(x$iso.n),]
    for(i in 2:nrow(x)){
      x$cor[i] = cor(r$intensity[x$id[1],],r$intensity[x$id[i],],m='sp')
      x$drt[i] = r$peaks[x$id[i],'rt']-r$peaks[x$id[1],'rt']
      x$dmz[i] = (r$peaks[x$id[i],'mz']-r$peaks[x$id[1],'mz'])/(i-1)*x$iso.charge[i] - 1.003355
    }
    x
  })
  
  iso = do.call(rbind,iso)
  rownames(iso) = iso$id
  iso = iso[rownames(r$peaks),]
  colnames(iso)[4:6] = paste0('iso.',colnames(iso)[4:6])
  r$peaks = cbind(r$peaks[,c('peak.id','rt','mz','isotopes')],iso[,1:6])
  r
}


generateFA = function(ns,ks,aw=atom.weights){
  r =NULL
  for(n in ns)
    for(k in ks){
      h = 2*n-2*k
      r = rbind(r,data.frame(LM_ID=paste0('FA_',n,':',k),
                             FORMULA=paste0('C',n,'H',h,'O2'),
                             n=n,k=k,
                             EXACT_MASS=aw['C','Mass']*n+aw['H','Mass']*h+aw['O','Mass']*2))
    }
  rownames(r) = r$LM_ID
  r
}

generateTGL = function(ns,ks,aw=atom.weights){
  r = NULL
  for(n in ns)
    for(k in ks)
      r = rbind(r,data.frame(n=n,k=k,FORMULA=paste0('C',3+n,'H',2*(n-k+1),'O6')))
    rownames(r)=r$FORMULA
    r$EXACT_MASS = sapply(r$FORMULA,calcExactMass,aws=aw)
    r$LM_ID = paste0('TAG_',r$n,':',r$k)
    r
}

generateTGLplus0 = function(ns,ks){
  r = NULL
  for(n in ns)
    for(k in ks){
      f = paste0('C',3+n,'H',2*(n-k+1),'O7')
      r = rbind(r,data.frame(LM_ID = paste0('TAG+O_',n,':',k),
                             EXACT_MASS = calcExactMass(f),
                             FORMULA=f,
                             sub.class='TAG+O',
                             n=n,k=k))
    }
    rownames(r)=r$LM_ID
    r
}

getPeaks = function(p,rt,mz)p[p$rt>rt[1] & p$rt<rt[2] & p$mz>mz[1] & p$mz<mz[2],]

cleanNetsByDirection = function(p){
	cleanNet = function(p,s){
		check = function(v) length(v) == 0 || sum(v) > 0 
		f = !is.na(p$start) & p$start == s
		for(i in which(f)){
			ok = TRUE
			ok = ok && check(p[f & p$k==p$k[i] & p$n == p$n[i]+1,'rt'] - p$rt[i] > 0)
			ok = ok && check(p[f & p$k==p$k[i] & p$n == p$n[i]-1,'rt'] - p$rt[i] < 0)
			ok = ok && check(p[f & p$n==p$n[i] & p$k == p$k[i]+1,'rt'] - p$rt[i] > 0)
			ok = ok && check(p[f & p$n==p$n[i] & p$k == p$k[i]-1,'rt'] - p$rt[i] < 0)
			if(!ok){
				p$start[i] = p$n[i] = p$k[i] = NA
				f[i] = FALSE
			}
		}
		p
	}
	#for group of peaks that have same net coordinates leave one these that have correct directions to all existing neighbors (at least on of them)
	gr = unique(na.omit(p$start))
	for(g in gr){
		repeat{
			f = is.na(p$start)
			p = cleanNet(p,g)
			if(all(f == is.na(p$start))) break
		}
	}
	p
}



annotateByNets = function(p,ann,ppm=20,plot=TRUE,clean.net.zs.thr=3,min.net.size=0){
  net.starts = sort(table(p$start),decreasing = T)
  net.starts = net.starts[net.starts>=min.net.size]
  net.no = 1
  p$net.no = NA
  p$pass = FALSE
  p$sub.class = NA
  for(ns in as.numeric(names(net.starts))){
    f1 = !is.na(p$start) & p$start == ns
    p$net.no[f1] = net.no
    if(sum(f1)>=0)
      cn = cleanNet(p[f1,],zs.thr = clean.net.zs.thr)
    else
      cn = p[f1,]
    
    f2 = f1 & rownames(p) %in% rownames(cn)
    m = lm(p$rt[f2] ~ (p$n[f2] + I(p$n[f2]^2)+ I(p$n[f2]^3))*(p$k[f2]+I(p$k[f2]^2)+I(p$k[f2]^3)))
    prt = predict(m,newdata = list(p=p,f2=f1))
   # if(!is.nan(summary(m)$sigma) && summary(m)$sigma>5/60) f2[] = FALSE
    p$pass[f2] = TRUE
    #find best annotation
    na = ann[ann$id %in% p$id[f2],]
    subs = sort(table(paste(na$sub.class,na$ion)),decreasing = T)
    if(nrow(na)>0){
      subs = subs[1:min(6,length(subs))]
      p$sub.class[f1] = names(subs)[1]
    }
    if(plot){
      plot(prt-p$rt[f1],col=ifelse(f2[f1],'green','red'),ylab='Residuals',main=paste0('Net ',net.no,' (',sum(f2),'/',sum(f1),')'))
      plot(p$rt[f1],p$mz[f1],t='n',xlab='RT',ylab='m/z')
      points(p$rt[f1 & !f2],p$mz[f1 & !f2],col='red')
      pch=1
      if(nrow(na)>0)
        pch=ifelse(p$id[f2] %in% na$id[paste(na$sub.class,na$ion)==names(subs)[1]],19,1)
      plotNet(p[f2,],FALSE,pch=pch)
    }
    #look on ann
    if(nrow(na)==0){
      plot.new()
      text(0.5,0.5,'No one peak\nhas annotation',cex=3)
    }else{
      ppms = split(na$ppm,paste(na$sub.class,na$ion))[names(subs)]
      boxplot(ppms,at=1:length(ppms),xaxt='n',notch = T,ylab='ppm',ylim=c(0,ppm),main=names(subs)[1])
      axis(1,1:length(ppms),paste0(names(subs),'\n',round(subs/sum(f2)*100),'%'),las=3)
      lines(as.numeric(subs)/sum(f2)*ppm,col='red',lwd=3,t='b')
      axis(4,0:5/5*ppm,0:5/5*100)
      mtext('% of annotated',4,1.3)
    }
    net.no = net.no + 1;
  }
  p
}

lookForNets = function(res,steps=c(n=14.01565,k=2.01565),step.max.cnt=c(n=2,k=1),ppm=20,rt.win = c(0,0.5)){
  # 0.5 of min is too small, there are nets with longer step... should I always use closest (by RT) node?
  # should I use threshold for minimal dRT?
  ord = order((1:nrow(res))[order(res$rt)])
  
  ppm = ppm*1e-6
  res$start=NA
  res = res[order(res$rt),]
  rt = res$rt
  mz = res$mz
  for(n in names(steps)) res[[n]] = NA
  for(i in 1:(nrow(res)-1)){
    #cat('\r',i)
    for(s in names(steps)){
      for(sc in 1:step.max.cnt[s]){
        #find RT window
        for(j in (i+1):nrow(res)) 
          if(rt[j]-rt[i] > rt.win[2]*sc){
            j = j - 1
            break
          }
      	
      	if(j > i)
      		minrt.i = i + which((rt[(i+1):j]-rt[i]) >= rt.win[1])[1]
      	else
      		minrt.i = NA

        if(!is.na(minrt.i) & j>=minrt.i){
          tmz = mz[minrt.i:j]
          ppms = abs(mz[i] + steps[s]*sc - tmz)/tmz
          mininx.=which.min(ppms)
          if(ppms[mininx.] < ppm){
            mininx = mininx. + minrt.i - 1
            if(is.na(res$start[i])){
              res$start[i] = i
              for(ss in names(steps)) res[[ss]][i] = 0
            }
            if(is.na(res$start[mininx])){
              res$start[mininx] = res$start[i]
              for(ss in names(steps)) res[[ss]][mininx] = res[[ss]][i]
              res[[s]][mininx] = res[[s]][i] + sc
            }else{
              #merge nets
              s1 = res$start[i]
              s2 = res$start[mininx]
              if(mz[s2] < mz[s1]){
                oldf = !is.na(res$start) & res$start==s1
                for(ss in names(steps))
                  res[[ss]][oldf] = res[[ss]][oldf]+res[[ss]][mininx]-res[[ss]][i]
                res[[s]][oldf] = res[[s]][oldf] - sc
                res$start[oldf] = s2
              }else{
                oldf = !is.na(res$start) & res$start==s2
                for(ss in names(steps))
                  res[[ss]][oldf] = res[[ss]][oldf]-res[[ss]][mininx]+res[[ss]][i]
                res[[s]][oldf] = res[[s]][oldf] + sc
                res$start[oldf] = s1
              }
            }
            break
          }
        }
      }
    }
  }
  res[ord,]
}

cleanNet = function(d,zs.thr,abs.thr=Inf){
  d$rownames = rownames(d)
  m = lm(d$rt ~ (d$n + I(d$n^2)+ I(d$n^3))*(d$k+I(d$k^2)+I(d$k^3)))
  r = residuals(m)
  d$r = abs((r-mean(r))/sd(r))
  f = d$r>zs.thr | abs(r) > abs.thr
  if(sum(f,na.rm=T)>0){
    d = d[!f,]
    d = cleanNet(d,zs.thr)
  }else{
    d=do.call(rbind,lapply(split(d,paste(d$n,d$k)),function(x){x[order(x$r)[1],]}))
    rownames(d) = d$rownames
    d$rownames = NULL
    d$r = NULL
  }
  d
}

generateNet = function(start.mass,ns,ks,base.name){
  r = NULL
  for(n in ns)
    for(k in ks)
      r = rbind(r,data.frame(n=n,k=k,FORMULA=''))
    r$LM_ID = paste0(base.name,r$n,':',r$k)
    rownames(r)=r$LM_ID
    r$EXACT_MASS = start.mass + r$n*calcExactMass('CH2') - r$k*calcExactMass('H2')
    r
}

plotNet = function(d,new=T,col=NULL,cex=NULL,col.lines='black',xlab='RT (min)',ylab='m/z',peak.labels=NULL,peak.labels.cex=1,mz.shift=0,...){
  d$mz = d$mz-mz.shift
  if(is.null(d$n))
    d$n = d$chain.len
  if(is.null(d$k))
    d$k = d$double.bonds
  if(is.null(col)){
    col = d$k-min(d$k) + 1
    #print(col)
    col = rainbow(max(col))[col]
  }
  if(is.null(cex)){
    if(max(d$n)!=min(d$n))
      cex = (d$n-min(d$n))/(max(d$n)-min(d$n))*1.8 + 0.8
    else
      cex = 1
  }
  if(new)
    plot(d$rt,d$mz,xlab=xlab,ylab=ylab,t='n',col=col,cex=cex,...)
  points(d$rt,d$mz,col=col,cex=cex,...)
  for(i in 1:nrow(d)){
    n = d$n[i]
    k = d$k[i]
    js = which(d$n == n-1 & d$k == k)
    if(length(js) == 0) js = which(d$n == n-2 & d$k == k)
    js = c(js,which(d$n == n   & d$k == k-1))
    segments(d$rt[js],d$mz[js],rep(d$rt[i],length(js)),rep(d$mz[i],length(js)),col=col.lines,...)
  }
  if(!is.null(peak.labels))
  	text(d$rt,d$mz,peak.labels,adj=c(0.2,-1),col=col,cex=peak.labels.cex)
}

getNegAdducts = function(){
  c('-H'=0,C2H4O2=calcExactMass('C2H4O2'),CH2O2=calcExactMass('CH2O2')) - calcExactMass('H') + 0.00054858
}

annotateByMass = function(mz.rt,db,ppm=100,delta=0.1,ions=c(H=1.007276,NH4=18.033823,Na=22.989218,K=38.96315,aNH4=calcExactMass('NH4')+calcExactMass('C2H3N')-0.00054858)){
  db = db[order(db$EXACT_MASS),]
  mz.rt = mz.rt[order(mz.rt$mz),]
  r=do.call(rbind,lapply(names(ions),function(i){annotateByMass.(mz.rt,db,ions[i],ppm,delta)}))
  r[order(r$mz,r$ppm,r$LM_ID),]
}

annotateByMass. = function(mz.rt,db,adduct,ppm=100,delta=0.1){
  r = new.env()
  j = 0;
  for(i in 1:nrow(mz.rt)){
    ##cat('\r',names(adduct),i,nrow(mz.rt),'     ')
    if(j == nrow(db)) break
    mo = mz.rt$mz[i]
    for(k in (j+1):nrow(db)){
      d = db$EXACT_MASS[k]+adduct-mo
      p = d/(db$EXACT_MASS[k]+adduct)*1e6
      if(ifelse(!is.null(ppm),abs(p)<=ppm,abs(d)<=delta)){
        if (dim(db)[2] == 5){
          db.colnames <- c('LM_ID','EXACT_MASS','FORMULA')
        }
        else{
          db.colnames <- c('LM_ID','EXACT_MASS','FORMULA', 'SYSTEMATIC_NAME', 'ABBREV')
        }
        r[[as.character(length(r))]] = cbind(mz.rt[i,], ion=names(adduct), db[k, db.colnames], ppm=abs(p), delta=abs(d), ppmd = p)
      }else{
        if(d < 0){
          j = k
          if(j == nrow(db)) break
        }
        else
          break
      }
    }
  }
  do.call(rbind,as.list(r))
}

calcExactMasses = function(fs,ret.formula=FALSE){
  sapply(fs,calcExactMass,ret.formula=ret.formula)
}

calcExactMass = function(f,aws=atom.weights,ret.formula=FALSE){
  if(is.na(f))
    return(NA)
  p = gregexpr('[A-Z][a-z]?\\d*',f,perl = TRUE)[[1]]
  l = attr(p,"match.length")
  w = 0
  
  for(i in 1:length(p)){
    e. = substr(f,p[i],p[i]+l[i]-1)
    e = gsub('\\d+','',e.)
    n = gsub('[A-Za-z]','',e.)
    if(n=='')
      n = 1
    else
      n=as.numeric(n)
    if(ret.formula)
      w[e] = n
    else
      w = w + n*aws[e,'Mass']
  }
  if(ret.formula)
    w = w[-1]
  w
}

mapPeaks2Peaks = function(p1,p2,max.dRT=2,max.dMZ=0.1){
  r = data.frame(inx1 = 1:nrow(p1))
  r$inx2 = r$dMZ = r$dRT = NA
  p2$inx=1:nrow(p2)
  for(i in 1:nrow(p1)){
    t = p2[abs(p2$rt-p1$rt[i]) <=max.dRT &  abs(p2$mz-p1$mz[i]) <=max.dMZ,]
    if(nrow(t)>0){
      t = t[order(abs(t$mz-p1$mz[i]))[1],]
      r$inx2[i] = t$inx
      r$dMZ[i] = t$mz-p1$mz[i]
      r$dRT[i] = t$rt-p1$rt[i]
    }
  }
  r
}

# MS2 ####
loadMS2 = function(f){
  require(mzR)
  x = openMSfile(f)
  runInfo(x)
  p=peaks(x)
  hs = header(x)
  table(hs$msLevel)
  hfs=c('msLevel','peaksCount','retentionTime','precursorMZ','precursorIntensity','basePeakMZ')
  hs. = hs[,hfs]
  close(x)
  list(peaks=p,header=hs.,header.full=hs)
}

getMaxMS2Peaks = function(s,n=5,dmz=0.05){
  f = s$header$msLevel==2
  h = s$header[f,]
  p = s$peaks[f]
  mz = int = matrix(NA,ncol=n,nrow=sum(f))
  for(i in 1:length(p)){
    t = p[[i]]
    for(k in 1:n){
      inx = which.max(t[,2])
      mz[i,k] = t[inx,1]
      int[i,k] = t[inx,2]
      t[abs(t[,1]-t[inx,1])<dmz,2] = 0
    }
  }
  list(prec=h,frag.mz=mz,frag.int=int)
}
