#include<iostream> //bcpt > kb.ppm
#include<cmath>
#define H return
using namespace std;typedef int N;typedef float F;N SPP=99;
F PI=3.1416;N WT=0;F E=1e-2;F FM=1e38;F Rg(){H(F)drand48();
}struct V{F x,y,z;V(F v):x(v),y(v),z(v){}V():x(0),y(0),z(0)
{}V(F x,F y,F z):x(x),y(y),z(z){}V operator-()const{H V(-x,
-y,-z);}F L(){H sqrt(*this%*this);}V operator+(const V&v)
const{H V(x+v.x,y+v.y,z+v.z);}V operator/(F s)const{H*this*
(1/s);}V operator*(const V&v)const{H V(x*v.x,y*v.y,z*v.z);}
V operator~(){H*this/L();}V operator*(F s)const{H V(x*s,y*s
,z*s);}V operator-(const V&v)const{H*this+(-v);}F operator%
(V&v){H v.x *x+v.y*y+v.z*z;}};V SS(){F z=1-Rg()*2,r=sqrt(
max(0.f,1-z*z)),phi=2*PI*Rg();H V(r*cos(phi),r*sin(phi),z);
}struct R{V o,d;F n=E,f=FM;};struct I{V p,n,k;F t=FM;N m=0;
};struct L{V p,i;};struct S{V c,k;F rd;N m;N I(R&r,I&i){V u
=r.o-c;F b=2*(r.d%u);F d=b*b-4*(u%u-rd*rd);if(d>0){d=sqrt(d
);F t=(-b-d)*.5f;if(t<=r.n)t=(-b+d)*.5f;if(t>r.n&&t<r.f){i.
p=r.o+r.d*t;i.n=~(i.p-c);i.t=t;i.k=k;i.m=m;H 1;}}H 0;}};L
lt[]={{V(-4,2,-2),V(200)},{V(0,2,0),V(110,100,80)}};S ob[59
]={{V(-505,0,0),V(1,.5,.5),500,1},{V(0,505,0),V(1),500,1},{
V(0,-5005,0),V(1),5000},{V(505,0,0),V(.5,1,.5),500,1},{V(0,
0,505),V(1),500,1},{V(2,2,-2),V(.9,.7,.6),1,2},{V(2,-4,-2),
V(.9,.9,.9),1,3},{-V(2,3,2),V(.9,.7,.6),2,2},{V(2,0,-2),V(
.9,.9,.9),1,3},{V(2,-2,-2),V(.6,.7,.9),1,2}};N In(R&r,I&i){
I z;N f=0;for(S&s:ob){if(s.I(r,z)){if(z.t<i.t)i=z;f=1;}}H f
;}N h[]={20,474,1498,3546,7700,3546,1498,474,20};N It(V&a,V
&b){I i;R r;r.d=b-a+SS();r.f=r.d.L()-E;r.n=E;r.o=a;r.d=~r.d
;for(S&s:ob)if(s.I(r,i))H 1;H 0;}V Fe(I&i){F k=(!i.m&&(N)(i
.p.x-99)%2!=(N)(i.p.z-99)%2)?.5f:1;H i.k*((i.m>1)?k:k/PI);}
V Fr(I&i,V&t,V&u){if(i.m<2){u=((u=SS())%i.n<0)?-u:u;}else
if(i.m<3){u=t-(i.n*(t%i.n*2));}else if(i.m<4){F r=2;V n=i.n
;F c=n%t;if (c>0){n=-n;c=-c;}else r=1/r;u=t*r-n*(r*c+sqrt(1
-r*r*(1-c*c)));}H Fe(i);}V Sh(I&i){V Le;for(L&l:lt)if(!It(i
.p,l.p)){V d=l.p-i.p;F a=d%d;d=~d;Le=Le+l.i*Fe(i)*max(0.f,d
%i.n)/a;}H Le;}V Tr(R&r,N b){I i;V le;V a(1);for(N bb=0;bb<
b;++bb){if(!In(r,i))break;V f=Fr(i,r.d,r.d);if(i.m<2){le=le
+Sh(i)*a;if(WT)break;a=a*abs(r.d%i.n);}r.o=i.p;r.n=E;r.f=i.
t=FM;a=a*f;}H le;}void Tm(F x){cout<<(uint)(x/(1 + x)*255)
<<' ';}N main(){N j=9;S a={{0,-4.8,-9},{.9,.8,.6},.2,2};for
(N i:h){i+=17e3;a.m^=1;for(N x=0;x<31;++x)if(i&2<<x){a.c.x=
(x-6)*.4f;ob[++j]=a;}a.c.z+=.4;}cout << "P3\n512 512 255\n"
;for(N y=-256;y<256;++y){for(N x=-224;x<288;++x){V l;for(N
s=0;s<SPP;++s){R r{-V(0,3,15),~(V(x+Rg(),y+Rg(),512)/512)};
l=l+Tr(r,20)/SPP;}Tm(l.x);Tm(l.y);Tm(l.z);}cout<<endl;}}