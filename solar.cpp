#include "iGraphics.h"
#include "gl.h"
#include<math.h>

double mx1,my1;
double pi=acos(-1)/180;

int line=1;
int space=1;

int sunx=500,suny=460;
float zoom=.70;

double merx,mery,merl=1,mt=90;
double venx,veny,venl=1,vt=25;
double earx,eary,earl=1,et=7;
    double moonx,moony,moonl,moont,ma,mb;
double marsx,marsy,marsl=1,marst=80;
    double m1x,m1y,m1l,m1t,m1a,m1b;
    double m2x,m2y,m2l,m2t,m2a,m2b;
double jx,jy,jl=1,jt=14;
    double sjx[5],sjy[5],sjl[5],sjt[5]={0,45,90,180,210},sja[5],sjb[5];
double sx,sy,sl=1,st=125;
    double ssx[3],ssy[3],ssl[3],sst[3]={210,300,100},ssa[3],ssb[3];
double ux,uy,ul=1,ut=180;
    double sux[2],suy[2],sul[2],sut[2]={0,210},sua[2],sub[2];
double nx,ny,nl=1,nt=210;
    double snx[2],sny[2],snl[2],snt[2]={45,180},sna[2],snb[2];

double a[8],b[8];
double cx[2]={-700,5100},cy[4]={800,-2300};

void iDraw(){
    iClear();

    //Background
    iShowBMP2(45,15,"stars.bmp",0);

    iShowBMP2(cx[0],cy[0],"tltobr.bmp",0);
    iShowBMP2(cx[1],cy[1],"brtotl.bmp",0);




    iSetColor(255,239,7);
    //SUN
    if(zoom<.7){
        iFilledCircle(sunx,suny,90*zoom);
    }
    else if(zoom>1.65)iShowBMP2(sunx-180,suny-200,"bigsun.bmp",0);
    else if(sunx-75<0 || suny-115<0){
        iFilledCircle(sunx,suny,90*zoom);
    }
    else iShowBMP2(sunx-75,suny-115,"sun.bmp",0);

    if(line){
        iSetColor(75,75,75);
        iEllipse(sunx+sqrt(a[0]*a[0]-b[0]*b[0]),suny,a[0],b[0]);//Mercury
        iEllipse(sunx+sqrt(a[1]*a[1]-b[1]*b[1]),suny,a[1],b[1]);//Venus
        iEllipse(sunx+sqrt(a[2]*a[2]-b[2]*b[2]),suny,a[2],b[2]);//Earth
             iSetColor(35,35,35);
             iEllipse(earx+sqrt(ma*ma-mb*mb),eary,ma,mb);//Moon
        iSetColor(75,75,75);
        iEllipse(sunx+sqrt(a[3]*a[3]-b[3]*b[3]),suny,a[3],b[3]);//Mars
            iSetColor(35,35,35);
            iEllipse(marsx+sqrt(m1a*m1a-m1b*m1b),marsy,m1a,m1b);//SubMars1
            iEllipse(marsx+sqrt(m2a*m2a-m2b*m2b),marsy,m2a,m2b);//SubMars2
        iSetColor(75,75,75);
        iEllipse(sunx+sqrt(a[4]*a[4]-b[4]*b[4]),suny,a[4],b[4]);//Jupiter
            iSetColor(35,35,35);
            iEllipse(jx+sqrt(sja[0]*sja[0]-sjb[0]*sjb[0]),jy,sja[0],sjb[0]);//SubJupiter1
            iEllipse(jx+sqrt(sja[1]*sja[1]-sjb[1]*sjb[1]),jy,sja[1],sjb[1]);//SubJupiter2
            iEllipse(jx+sqrt(sja[2]*sja[2]-sjb[2]*sjb[2]),jy,sja[2],sjb[2]);//SubJupiter3
            iEllipse(jx+sqrt(sja[3]*sja[3]-sjb[3]*sjb[3]),jy,sja[3],sjb[3]);//SubJupiter4
            iEllipse(jx+sqrt(sja[4]*sja[4]-sjb[4]*sjb[4]),jy,sja[4],sjb[4]);//SubJupiter5
        iSetColor(75,75,75);
        iEllipse(sunx+sqrt(a[5]*a[5]-b[5]*b[5]),suny,a[5],b[5]);//Saturn
            iSetColor(35,35,35);
            iEllipse(sx+sqrt(ssa[0]*ssa[0]-ssb[0]*ssb[0]),sy,ssa[0],ssb[0]);//SubSaturn1
            iEllipse(sx+sqrt(ssa[1]*ssa[1]-ssb[1]*ssb[1]),sy,ssa[1],ssb[1]);//SubSaturn2
            iEllipse(sx+sqrt(ssa[2]*ssa[2]-ssb[2]*ssb[2]),sy,ssa[2],ssb[2]);//SubSaturn3
        iSetColor(75,75,75);
        iEllipse(sunx+sqrt(a[6]*a[6]-b[6]*b[6]),suny,a[6],b[6]);//Uranus
            iSetColor(35,35,35);
            iEllipse(ux+sqrt(sua[0]*sua[0]-sub[0]*sub[0]),uy,sua[0],sub[0]);//SubUranus1
            iEllipse(ux+sqrt(sua[1]*sua[1]-sub[1]*sub[1]),uy,sua[1],sub[1]);//SubUranus2
        iSetColor(75,75,75);
        iEllipse(sunx+sqrt(a[7]*a[7]-b[7]*b[7]),suny,a[7],b[7]);//Neptune
            iSetColor(35,35,35);
            iEllipse(nx+sqrt(sna[0]*sna[0]-snb[0]*snb[0]),ny,sna[0],snb[0]);//SubNeptune1
            iEllipse(nx+sqrt(sna[1]*sna[1]-snb[1]*snb[1]),ny,sna[1],snb[1]);//SubNeptune2
    }
    //MERCURY
    iSetColor(255,255,204);
    if(zoom<.7)
        iFilledCircle(merx,mery,10*zoom);
    else if(zoom>1.65)iShowBMP2(merx-14,mery-14,"bigmercury.bmp",0);
    else iShowBMP2(merx-7,mery-7,"mercury.bmp",0);

    //VENUS
    iSetColor(255,178,102);
    if(zoom<.7)
        iFilledCircle(venx,veny,11*zoom);
    else if(zoom>1.65)iShowBMP2(venx-18,veny-18,"bigvenus.bmp",0);
    else iShowBMP2(venx-9,veny-9,"venus.bmp",0);

    //EARTH
    iSetColor(51,51,255);
    if(zoom<.7)
        iFilledCircle(earx,eary,15*zoom);
    else if(zoom>1.65)
        iShowBMP2(earx-26,eary-26,"bigearth.bmp",0);
    else
        iShowBMP2(earx-13,eary-13,"earth.bmp",0);
        //MOON
        iSetColor(250,245,255);
        if(zoom<.7)iFilledCircle(moonx,moony,7*zoom);
        else if(zoom>1.65)iShowBMP2(moonx-10,moony-10,"bigmoon.bmp",0);
        else iShowBMP2(moonx-5,moony-5,"moon.bmp",0);

    //MARS
    iSetColor(255,204,153);
    if(zoom<.7)
        iFilledCircle(marsx,marsy,13*zoom);
    else if(zoom>1.65)iShowBMP2(marsx-20,marsy-20,"bigmars.bmp",0);
    else iShowBMP2(marsx-10,marsy-10,"mars.bmp",0);
        //Sub-Mars
        iSetColor(254,169,43);
        if(zoom<.7)iFilledCircle(m1x,m1y,6*zoom);
        else if(zoom>1.65)iShowBMP2(m1x-7,m1y-7,"bigsmars1.bmp",0);
        else iShowBMP2(m1x-3.5,m1y-3.5,"smars1.bmp",0);
        iSetColor(24,143,255);
        if(zoom<.7)iFilledCircle(m2x,m2y,7*zoom);
        else if(zoom>1.65)iShowBMP2(m2x-8,m2y-8,"bigsmars2.bmp",0);
        else iShowBMP2(m2x-4,m2y-4,"smars2.bmp",0);

    //JUPITER
    iSetColor(255,178,102);
    if(zoom<.7)iFilledCircle(jx,jy,35*zoom);
    else if(zoom>1.65)iShowBMP2(jx-50,jy-50,"bigjupiter.bmp",0);
    else iShowBMP2(jx-25,jy-25,"jupiter.bmp",0);
        //Sb-Planet
        iSetColor(254,169,43);
        if(zoom<.7)iFilledCircle(sjx[0],sjy[0],7*zoom);
        else if(zoom>1.65)iShowBMP2(sjx[0]-11,sjy[0]-11,"bigsubjupiter1.bmp",0);
        else iShowBMP2(sjx[0]-5.5,sjy[0]-5.5,"subjupiter1.bmp",0);
        iSetColor(24,143,255);
        if(zoom<.7)iFilledCircle(sjx[1],sjy[1],4.5*zoom);
        else if(zoom>1.65)iShowBMP2(sjx[1]-6,sjy[1]-7,"bigsubjupiter2.bmp",0);
        else iShowBMP2(sjx[1]-3.5,sjy[1]-3.5,"subjupiter2.bmp",0);
        iSetColor(204,143,255);
        if(zoom<.7)iFilledCircle(sjx[2],sjy[2],5.5*zoom);
        else if(zoom>1.65)iShowBMP2(sjx[2]-6,sjy[2]-6,"bigsmars1.bmp",0);
        else iShowBMP2(sjx[2]-3,sjy[2]-3,"smars1.bmp",0);
        iSetColor(24,243,55);
        if(zoom<.7)iFilledCircle(sjx[3],sjy[3],7.5*zoom);
        else if(zoom>1.65)iShowBMP2(sjx[3]-3,sjy[3]-3,"bigsmars2.bmp",0);
        else iShowBMP2(sjx[3]-3,sjy[3]-3,"smars2.bmp",0);
        iSetColor(45,201,156);
        if(zoom<.7)iFilledCircle(sjx[4],sjy[4],9*zoom);
        else if(zoom>1.65)iShowBMP2(sjx[4]-3,sjy[4]-3,"bigsubsaturn1.bmp",0);
        else iShowBMP2(sjx[4]-3,sjy[4]-3,"subsaturn1.bmp",0);

    //SATURN
    iSetColor(255,229,204);
    if(zoom<.7){
        iFilledCircle(sx,sy,30*zoom);
        iSetColor(255,218,200);
        iEllipse(sx,sy,40*zoom,20*zoom);
        iEllipse(sx,sy,39*zoom,19*zoom);
        iEllipse(sx,sy,38*zoom,18*zoom);
        iEllipse(sx,sy,37*zoom,17*zoom);
    }
    else if(zoom>1.65){
        if(sx>=sunx+sqrt(a[5]*a[5]-b[5]*b[5])*2)
            iShowBMP2(sx-40,sy-40,"bigsatright.bmp",0);
        else if(sx>=sunx)
            iShowBMP2(sx-56,sy-30,"bigsatstr.bmp",0);
        else
            iShowBMP2(sx-40,sy-40,"bigsatleft.bmp",0);
    }
    else{
        if(sx>=sunx+sqrt(a[5]*a[5]-b[5]*b[5])*2)
            iShowBMP2(sx-20,sy-20,"satright2.bmp",0);
        else if(sx>=sunx)
            iShowBMP2(sx-28,sy-15,"satstr.bmp",0);
        else
            iShowBMP2(sx-20,sy-20,"satleft2.bmp",0);
    }
    //Sub-Saturn
        iSetColor(25,169,43);
        if(zoom<.7)iFilledCircle(ssx[0],ssy[0],7*zoom);
        else if(zoom>1.65)iShowBMP2(ssx[0]-14,ssy[0]-14,"bigsubsaturn1.bmp",0);
        else iShowBMP2(ssx[0]-7,ssy[0]-7,"subsaturn1.bmp",0);
        iSetColor(24,14,255);
        if(zoom<.7)iFilledCircle(ssx[1],ssy[1],11*zoom);
        else if(zoom>1.65)iShowBMP2(ssx[1]-22,ssy[1]-22,"bigsubsaturn2.bmp",0);
        else iShowBMP2(ssx[1]-11,ssy[1]-11,"subsaturn2.bmp",0);
        iSetColor(204,143,255);
        if(zoom<.7)iFilledCircle(ssx[2],ssy[2],9*zoom);
        else if(zoom>1.65)iShowBMP2(ssx[2]-18,ssy[2]-18,"bigsubsaturn3.bmp",0);
        else iShowBMP2(ssx[2]-9,ssy[2]-9,"subsaturn3.bmp",0);

    //URANUS
    iSetColor(51,51,255);
    if(zoom<.7)iFilledCircle(ux,uy,25*zoom);
    else if(zoom>1.65)iShowBMP2(ux-40,uy-40,"biguranus.bmp",0);
    else iShowBMP2(ux-20,uy-20,"uranus.bmp",0);
        //SubUranus
        iSetColor(254,169,43);
        if(zoom<.7)iFilledCircle(sux[0],suy[0],8.5*zoom);
        else if(zoom>1.65)iShowBMP2(sux[0]-5.5,suy[0]-5.5,"bigsubjupiter1.bmp",0);
        else iShowBMP2(sux[0]-5.5,suy[0]-5.5,"subjupiter1.bmp",0);
        iSetColor(24,143,255);
        if(zoom<.7)iFilledCircle(sux[1],suy[1],6.5*zoom);
        else if(zoom>1.65)iShowBMP2(sux[1]-3,suy[1]-3,"bigsubjupiter2.bmp",0);
        else iShowBMP2(sux[1]-3,suy[1]-3,"subjupiter2.bmp",0);

    //NEPTUNE
    iSetColor(51,153,255);
    if(zoom<.7)iFilledCircle(nx,ny,22.5*zoom);
    else if(zoom>1.65)iShowBMP2(nx-45,ny-45,"bigneptune.bmp",0);
    else iShowBMP2(nx-22.5,ny-22.5,"neptune.bmp",0);
    //SubNeptune
        iSetColor(254,19,43);
        if(zoom<.7)iFilledCircle(snx[0],sny[0],7.5*zoom);
        else if(zoom>1.65)iShowBMP2(snx[0]-5.5,sny[0]-5.5,"bigsubjupiter1.bmp",0);
        else iShowBMP2(snx[0]-5.5,sny[0]-5.5,"subjupiter1.bmp",0);
        iSetColor(240,143,255);
        if(zoom<.7)iFilledCircle(snx[1],sny[1],5.5*zoom);
        else if(zoom>1.65)iShowBMP2(snx[1]-3,sny[1]-3,"bigsubjupiter2.bmp",0);
        else iShowBMP2(snx[1]-3,sny[1]-3,"subjupiter2.bmp",0);

        iText(100,70,"Press space to see control instruction",GLUT_BITMAP_TIMES_ROMAN_24);

    if(!space){
        iSetColor(191,191,191);
        iFilledRectangle(450,0,600,900);
        iSetColor(rand()%255,rand()%255,rand()%255);
        iText(655,750,"Press p pause",GLUT_BITMAP_TIMES_ROMAN_24);
        iText(629,720,"Press r for resume",GLUT_BITMAP_TIMES_ROMAN_24);
        iText(615,690,"Press o for orbits path",GLUT_BITMAP_TIMES_ROMAN_24);
        iText(547,660,"Use navigation key to set position",GLUT_BITMAP_TIMES_ROMAN_24);
        iText(628.5,630,"Press + to zoom in",GLUT_BITMAP_TIMES_ROMAN_24);
        iText(621,600,"Press - to zoom out",GLUT_BITMAP_TIMES_ROMAN_24);
        iText(542,480,"Mouse can be used to change position",GLUT_BITMAP_TIMES_ROMAN_24);
        iText(715,450,"and",GLUT_BITMAP_TIMES_ROMAN_24);
        iText(650,420,"zoom in or out",GLUT_BITMAP_TIMES_ROMAN_24);
        iText(615,350,"Press space for go back",GLUT_BITMAP_TIMES_ROMAN_24);
    }
}
void Planet(){

    //r=(b[0]*b[0])/((a+cos(mt*pi)*sqrt(a[0]*a[0]-b[0]*b[0]));
    a[0]=280*zoom;//Mercury
    b[0]=220*zoom;
    merl=(b[0]*b[0])/(a[0]+cos(mt*pi)*sqrt((a[0]*a[0])-(b[0]*b[0])));
    merx=sunx+sqrt(a[0]*a[0]-b[0]*b[0])*2+merl*cos(mt*pi);
    mery=suny+merl*sin(mt*pi);

//Venus
    a[1]=320*zoom;//Venus
    b[1]=252*zoom;
    venl=(b[1]*b[1])/(a[1]+cos(vt*pi)*sqrt(a[1]*a[1]-b[1]*b[1]));
    venx=sunx+sqrt(a[1]*a[1]-b[1]*b[1])*2+venl*cos(vt*pi);
    veny=suny+venl*sin(vt*pi);

//Earth
    a[2]=360*zoom;//Earth
    b[2]=287.5*zoom;
    earl=pow(b[2],2)/(a[2]+cos(et*pi)*sqrt(a[2]*a[2]-b[2]*b[2]));
    earx=sunx+sqrt(a[2]*a[2]-b[2]*b[2])*2+earl*cos(et*pi);
    eary=suny+earl*sin(et*pi);

//Mars
    //r=b^2/((a+cos(theta)*(a^2-b^2)^(1/2))
    a[3]=400*zoom;//Mars
    b[3]=325*zoom;
    marsl=pow(b[3],2)/(a[3]+cos(marst*pi)*sqrt(a[3]*a[3]-b[3]*b[3]));
    marsx=sunx+sqrt(a[3]*a[3]-b[3]*b[3])*2+marsl*cos(marst*pi);
    marsy=suny+marsl*sin(marst*pi);

//Jupiter
    //r=b^2/((a+cos(theta)*(a^2-b^2)^(1/2))
    a[4]=460*zoom;//Jupiter
    b[4]=380*zoom;
    jl=pow(b[4],2)/(a[4]+cos(jt*pi)*sqrt(a[4]*a[4]-b[4]*b[4]));
    jx=sunx+sqrt(a[4]*a[4]-b[4]*b[4])*2+jl*cos(jt*pi);
    jy=suny+jl*sin(jt*pi);

//Saturn
    //r=b^2/((a+cos(theta)*(a^2-b^2)^(1/2))
    a[5]=520*zoom;//Saturn
    b[5]=435.5*zoom;
    sl=pow(b[5],2)/(a[5]+cos(st*pi)*sqrt(a[5]*a[5]-b[5]*b[5]));
    sx=sunx+sqrt(a[5]*a[5]-b[5]*b[5])*2+sl*cos(st*pi);
    sy=suny+sl*sin(st*pi);

//Uranus
    //r=b^2/((a+cos(theta)*(a^2-b^2)^(1/2))
    a[6]=590*zoom;//Uranus
    b[6]=501.5*zoom;
    ul=pow(b[6],2)/(a[6]+cos(ut*pi)*sqrt(a[6]*a[6]-b[6]*b[6]));
    ux=sunx+sqrt(a[6]*a[6]-b[6]*b[6])*2+ul*cos(ut*pi);
    uy=suny+ul*sin(ut*pi);

//Neptune
    //r=b^2/((a+cos(theta)*(a^2-b^2)^(1/2))
    a[7]=675*zoom;//Neptune
    b[7]=582.1875*zoom;
    nl=pow(b[7],2)/(a[7]+cos(nt*pi)*sqrt(a[7]*a[7]-b[7]*b[7]));
    nx=sunx+sqrt(a[7]*a[7]-b[7]*b[7])*2+nl*cos(nt*pi);
    ny=suny+nl*sin(nt*pi);

}

void angle(){
//Mercury
    mt+=1.5;
    if(mt>360){
        mt-=360;
    }
//Venus
    vt+=1.35;
    if(vt>360){
        vt-=360;
    }
//Earth
    et+=1;
    if(et>360){
        et-=360;
    }
//Mars
    marst+=.9;
    if(marst>360){
        marst-=360;
    }
//Jupiter
    jt+=.8;
    if(jt>360){
        jt-=360;
    }
//Saturn
    st+=.75;
    if(st>360){
        st-=360;
    }
//Uranus
    ut+=.75;
    if(ut>360){
        ut-=360;
    }
//Neptune
    nt+=.6;
    if(nt>360){
        nt-=360;
    }


//Moon
    moont+=9;
        if(moont>360){
            moont-=360;
        }
//Sub-Mars
         m1t+=7.5;
        if(m1t>360){
        m1t-=360;
        }
         m2t+=6;
        if(m2t>360){
        m2t-=360;
        }
//Sub-Jupiter
        sjt[0]+=6;
        if(sjt[0]>360){
            sjt[0]-=360;
        }
        sjt[1]+=5.5;
        if(sjt[1]>360){
        sjt[1]-=360;
        }
         sjt[2]+=6.2;
        if(sjt[2]>360){
        sjt[2]-=360;
        }
        sjt[3]+=5.5;
        if(sjt[3]>360){
        sjt[3]-=360;
        }
        sjt[4]+=4.75;
        if(sjt[4]>360){
        sjt[4]-=360;
        }
//Sub-Saturn
        sst[0]+=5.2;
        if(sst[0]>360){
            sst[0]-=360;
        }
        sst[1]+=5.9;
        if(sst[1]>360){
        sst[1]-=360;
        }
         sst[2]+=5.65;
        if(sst[2]>360){
        sst[2]-=360;
        }
//Sub-Uranus
        sut[0]+=5.4;
        if(sut[0]>360){
            sut[0]-=360;
        }
        sut[1]+=6.75;
        if(sut[1]>360){
        sut[1]-=360;
        }
//Sub-Neptune
        snt[0]+=6.5;
        if(snt[0]>360){
            snt[0]-=360;
        }
        snt[1]+=7.5;
        if(snt[1]>360){
        snt[1]-=360;
        }

}

void subplanet(){

    //Earth
        ma=40*zoom;
        mb=35*zoom;
        moonl=pow(mb,2)/(ma+cos(pi*moont)*sqrt(ma*ma-mb*mb));
        moonx=earx+sqrt(ma*ma-mb*mb)*2+moonl*cos(pi*moont);
        moony=eary+moonl*sin(pi*moont);

    //Mars
        m1a=30*zoom;
        m1b=25*zoom;
        m1l=pow(m1b,2)/(m1a+cos(pi*m1t)*sqrt(m1a*m1a-m1b*m1b));
        m1x=marsx+sqrt(m1a*m1a-m1b*m1b)*2+m1l*cos(pi*m1t);
        m1y=marsy+m1l*sin(pi*m1t);
        m2a=35*zoom;
        m2b=30*zoom;
        m2l=pow(m2b,2)/(m2a+cos(pi*m2t)*sqrt(m2a*m2a-m2b*m2b));
        m2x=marsx+sqrt(m2a*m2a-m2b*m2b)*2+m2l*cos(pi*m2t);
        m2y=marsy+m2l*sin(pi*m2t);

    //Jupiter
        sja[0]=60*zoom;
        sjb[0]=55*zoom;
        sjl[0]=pow(sjb[0],2)/(sja[0]+cos(pi*sjt[0])*sqrt(sja[0]*sja[0]-sjb[0]*sjb[0]));
        sjx[0]=jx+sqrt(sja[0]*sja[0]-sjb[0]*sjb[0])*2+sjl[0]*cos(pi*sjt[0]);
        sjy[0]=jy+sjl[0]*sin(pi*sjt[0]);
        sja[1]=70*zoom;
        sjb[1]=62*zoom;
        sjl[1]=pow(sjb[1],2)/(sja[1]+cos(pi*sjt[1])*sqrt(sja[1]*sja[1]-sjb[1]*sjb[1]));
        sjx[1]=jx+sqrt(sja[1]*sja[1]-sjb[1]*sjb[1])*2+sjl[1]*cos(pi*sjt[1]);
        sjy[1]=jy+sjl[1]*sin(pi*sjt[1]);
        sja[2]=72*zoom;
        sjb[2]=65*zoom;
        sjl[2]=pow(sjb[2],2)/(sja[2]+cos(pi*sjt[2])*sqrt(sja[2]*sja[2]-sjb[2]*sjb[2]));
        sjx[2]=jx+sqrt(sja[2]*sja[2]-sjb[2]*sjb[2])*2+sjl[2]*cos(pi*sjt[2]);
        sjy[2]=jy+sjl[2]*sin(pi*sjt[2]);
        sja[3]=75*zoom;
        sjb[3]=65*zoom;
        sjl[3]=pow(sjb[3],2)/(sja[3]+cos(pi*sjt[3])*sqrt(sja[3]*sja[3]-sjb[3]*sjb[3]));
        sjx[3]=jx+sqrt(sja[3]*sja[3]-sjb[3]*sjb[3])*2+sjl[3]*cos(pi*sjt[3]);
        sjy[3]=jy+sjl[3]*sin(pi*sjt[3]);
        sja[4]=68*zoom;
        sjb[4]=62*zoom;
        sjl[4]=pow(sjb[4],2)/(sja[4]+cos(pi*sjt[4])*sqrt(sja[4]*sja[4]-sjb[4]*sjb[4]));
        sjx[4]=jx+sqrt(sja[4]*sja[4]-sjb[4]*sjb[4])*2+sjl[4]*cos(pi*sjt[4]);
        sjy[4]=jy+sjl[4]*sin(pi*sjt[4]);

    //Saturn
        ssa[0]=60*zoom;
        ssb[0]=55*zoom;
        ssl[0]=pow(ssb[0],2)/(ssa[0]+cos(pi*sst[0])*sqrt(ssa[0]*ssa[0]-ssb[0]*ssb[0]));
        ssx[0]=sx+sqrt(ssa[0]*ssa[0]-ssb[0]*ssb[0])*2+ssl[0]*cos(pi*sst[0]);
        ssy[0]=sy+ssl[0]*sin(pi*sst[0]);
        ssa[1]=70*zoom;
        ssb[1]=62*zoom;
        ssl[1]=pow(ssb[1],2)/(ssa[1]+cos(pi*sst[1])*sqrt(ssa[1]*ssa[1]-ssb[1]*ssb[1]));
        ssx[1]=sx+sqrt(ssa[1]*ssa[1]-ssb[1]*ssb[1])*2+ssl[1]*cos(pi*sst[1]);
        ssy[1]=sy+ssl[1]*sin(pi*sst[1]);
        ssa[2]=72*zoom;
        ssb[2]=65*zoom;
        ssl[2]=pow(ssb[2],2)/(ssa[2]+cos(pi*sst[2])*sqrt(ssa[2]*ssa[2]-ssb[2]*ssb[2]));
        ssx[2]=sx+sqrt(ssa[2]*ssa[2]-ssb[2]*ssb[2])*2+ssl[2]*cos(pi*sst[2]);
        ssy[2]=sy+ssl[2]*sin(pi*sst[2]);
    //Uranus
        sua[0]=60*zoom;
        sub[0]=55*zoom;
        sul[0]=pow(sub[0],2)/(sua[0]+cos(pi*sut[0])*sqrt(sua[0]*sua[0]-sub[0]*sub[0]));
        sux[0]=ux+sqrt(sua[0]*sua[0]-sub[0]*sub[0])*2+sul[0]*cos(pi*sut[0]);
        suy[0]=uy+sul[0]*sin(pi*sut[0]);
        sua[1]=70*zoom;
        sub[1]=62*zoom;
        sul[1]=pow(sub[1],2)/(sua[1]+cos(pi*sut[1])*sqrt(sua[1]*sua[1]-sub[1]*sub[1]));
        sux[1]=ux+sqrt(sua[1]*sua[1]-sub[1]*sub[1])*2+sul[1]*cos(pi*sut[1]);
        suy[1]=uy+sul[1]*sin(pi*sut[1]);
    //Neptune
        sna[0]=75*zoom;
        snb[0]=60*zoom;
        snl[0]=pow(snb[0],2)/(sna[0]+cos(pi*snt[0])*sqrt(sna[0]*sna[0]-snb[0]*snb[0]));
        snx[0]=nx+sqrt(sna[0]*sna[0]-snb[0]*snb[0])*2+snl[0]*cos(pi*snt[0]);
        sny[0]=ny+snl[0]*sin(pi*snt[0]);
        sna[1]=70*zoom;
        snb[1]=62*zoom;
        snl[1]=pow(snb[1],2)/(sna[1]+cos(pi*snt[1])*sqrt(sna[1]*sna[1]-snb[1]*snb[1]));
        snx[1]=nx+sqrt(sna[1]*sna[1]-snb[1]*snb[1])*2+snl[1]*cos(pi*snt[1]);
        sny[1]=ny+snl[1]*sin(pi*snt[1]);
}

void comet(){
    cx[0]+=100;
    cy[0]-=70;
    if(cy[0]<-10000){
        cx[0]=-4100;
        cy[0]=3220;
    }
    cx[1]-=100;
    cy[1]+=70;
    if(cy[1]>10860){
        cy[1]=-2300;
        cx[1]=5100;
    }

}

void iKeyboard(unsigned char key){
    if(key=='o'){
        line =!line;
    }
    if(key=='p'){
        iPauseTimer(3);
        iPauseTimer(1);

    }
    if(key=='r'){
        iResumeTimer(3);
        iResumeTimer(1);
    }
    if(key == '+'){
        if(zoom<2.2)
        zoom +=.05;
    }
    if(key=='-'){
        if(zoom>.45)
        zoom-=.05;
    }
    if(key==' '){
        space=!space;
    }
}

void iSpecialKeyboard(unsigned char key){
    if(key == GLUT_KEY_RIGHT){
        sunx -= 15;
    }
        if(key == GLUT_KEY_LEFT){
        sunx += 15;
    }
        if(key == GLUT_KEY_UP){
        suny -= 15;
    }
        if(key == GLUT_KEY_DOWN){
        suny += 15;
    }

}

void iMouse(int mouse,int state,int mx,int my)
{
    if(mouse == GLUT_RIGHT_BUTTON && state == GLUT_UP){
        if(zoom>.45)
        zoom-=.05;
    }
    if(mouse == GLUT_LEFT_BUTTON && state == GLUT_UP){
        if(zoom<2.2)
        zoom+=.05;
    }
}

void iMouseMove(int mx,int my){
    sunx=mx;
    suny=my;
}

int main(){
    iSetTimer(35,Planet);

    iSetTimer(35,angle);

    iSetTimer(25,subplanet);

    iSetTimer(100,comet);

    iInitialize(1600,900,"Solar System");
}
