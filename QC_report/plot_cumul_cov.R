# =======
#   License
# =======
#   This code is released under the GNU General Public License 3.0. A copy
# of this license is in the LICENSE.txt file.
# copyright Irina Krier 2015
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

plot_cumul_cov=function(cov,name,plot=FALSE,mincov=0)
{
  cov=cov[cov$index==name,]
  if(min(cov$V5)>mincov)
  {
    cov=rbind(cov[1,],cov)
    cov[1,c(5,6,8)]=c(0,0,0)
  }
  cov_cumul <- 1-cumsum(cov[,8])
  cols <- brewer.pal(3, "Dark2")
  col=cols[1]
  if(plot)
  {
    par(oma=c(1,2.8,2,0))
    plot(cov[2:nrow(cov), 5], cov_cumul[1:(length(cov_cumul)-1)], type='l',lwd=3,col=col, 
         ylim=c(0,1.0),xlim=c(0,60000),bty="n",axes=F,xlab="",ylab="")
    axis(1,cex.axis=0.7)
    axis(2,las=2,cex.axis=0.7)
    mtext(3,1,text=paste("Target Region Coverage :",format(100-100*cov$V8[cov$V5<=mincov],digits=2),"%"))
    mtext("Depth",1,2,cex=0.7)
    mtext("Fraction of capture target bases \u2265 depth",2,2,cex=0.7)
  }
  return(100-100*sum(cov$V8[cov$V5<=mincov]))
}