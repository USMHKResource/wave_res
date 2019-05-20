import wave_res as wr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
import shapefile as shp
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def setExtent(ex,region,ax,ret=False):
	
	if not ret:
		if region is 'ak':
			ex = [ex[0]-2.5,ex[1]+2,ex[2]-5,ex[3]]
			ax.set_extent(ex,crs=ccrs.PlateCarree())
		elif region is 'wc':
			ex = [ex[0]-5,ex[1],ex[2]-1,ex[3]]
			ax.set_extent(ex,crs=ccrs.PlateCarree())

		else:
			ax.set_extent(ex,crs=ccrs.PlateCarree())
	else:
		if region is 'ak':
			ex = [ex[0]-2.5,ex[1]+2,ex[2]-5,ex[3]]
		elif region is 'wc':
			ex = [ex[0]-5,ex[1],ex[2]-1,ex[3]]
		else:
			ex = ex
		
		return ex
	
		

def vectors(region,extent,ax,localArrows=False):
	ar = np.array
	if region is 'ak':
		lon, lat = ar([np.mean(extent[0:2])]), ar([np.mean([extent[-2:]])])
		
		# Remote Arrows
		lonShift,latShift = ([-29,-15,-5,5,15],
							 [-12,-13,-12,-10,-7])
		uf,vf = [10,10,10,10,30],[18,35,40,75,90]
		for i,shift in enumerate(zip(lonShift,latShift)):
			u, v = ar([uf[i]]),ar([vf[i]])
			ax.quiver(lon+shift[0],lat+shift[1],u,v,transform=ccrs.PlateCarree(),
					 minlength=10,width=0.005,headlength=5,pivot='middle')
		
		# Local Arrows
		if localArrows:
			lonShift,latShift = ([-20,-5,-5,15],
								 [-5,5,-6,-1])
			uf,vf = [-160,20,-40,50],[18,35,-40,75]
			for i,shift in enumerate(zip(lonShift,latShift)):
				u, v = ar([uf[i]]),ar([vf[i]])
				ax.quiver(lon+shift[0],lat+shift[1],u,v,transform=ccrs.PlateCarree(),
						 minlength=10,width=0.005,headlength=5,pivot='middle')
	
	elif region is 'at':
		lon, lat = ar([np.mean(extent[0:2])]), ar([np.mean([extent[-2:]])])
		
		# Remote Arrows
		lonShift,latShift = ([-12,-8,10,12,15.5],
							 [-10,-10.5,-5,-1,2.5])
		uf,vf = [10,10,-60,-70,-80],[18,35,40,75,90]
		for i,shift in enumerate(zip(lonShift,latShift)):
			u, v = ar([uf[i]]),ar([vf[i]])
			ax.quiver(lon+shift[0],lat+shift[1],u,v,transform=ccrs.PlateCarree(),
					 minlength=10,width=0.005,headlength=5,pivot='middle')
		
		# Local Arrows
		if localArrows:
			lonShift,latShift = [5,7,10.5], [-2,1,4.5]
			uf,vf = [-160,20,-40,50],[18,35,-40,75]
			for i,shift in enumerate(zip(lonShift,latShift)):
				u, v = ar([uf[i]]),ar([vf[i]])
				ax.quiver(lon+shift[0],lat+shift[1],u,v,transform=ccrs.PlateCarree(),
						 minlength=10,width=0.005,headlength=5,pivot='middle')
	elif region is 'prusvi':
		lon, lat = ar([np.mean(extent[0:2])]), ar([np.mean([extent[-2:]])])
		
		# Remote Arrows
		lonShift,latShift = [-3,-3,1,3,3],[-3,3,-3,-2.5,3]
		uf,vf = [10,10,-60,-70,-80],[18,-25,40,75,-60]
		for i,shift in enumerate(zip(lonShift,latShift)):
			u, v = ar([uf[i]]),ar([vf[i]])
			ax.quiver(lon+shift[0],lat+shift[1],u,v,transform=ccrs.PlateCarree(),
					 minlength=10,width=0.005,headlength=5,pivot='middle')
		
		# Local Arrows
		if localArrows:
			lonShift,latShift = [-0.5,0.75,0], [-2,-1.5,2]
			uf,vf = [-160,20,-40,50],[18,35,-40,75]
			for i,shift in enumerate(zip(lonShift,latShift)):
				u, v = ar([uf[i]]),ar([vf[i]])
				ax.quiver(lon+shift[0],lat+shift[1],u,v,transform=ccrs.PlateCarree(),
						 minlength=10,width=0.005,headlength=5,pivot='middle')
	elif region is 'wc':
		lon, lat = ar([np.mean(extent[0:2])]), ar([np.mean([extent[-2:]])])
		
		# Remote Arrows
		lonShift,latShift = [-10,-10,-9,-7,-5],[5,2,-1,-4,-7]
		uf,vf = [10,10,10,10,10],[1,3,8,9,12]
		for i,shift in enumerate(zip(lonShift,latShift)):
			u, v = ar([uf[i]]),ar([vf[i]])
			ax.quiver(lon+shift[0],lat+shift[1],u,v,transform=ccrs.PlateCarree(),
					 minlength=10,width=0.005,headlength=5,pivot='middle')
		
		# Local Arrows
		if localArrows:
			lonShift,latShift = [-3,-3,-2,0],[2,0,-3,-6]
			uf,vf = [-160,20,-40,50],[18,35,-40,75]
			for i,shift in enumerate(zip(lonShift,latShift)):
				u, v = ar([uf[i]]),ar([vf[i]])
				ax.quiver(lon+shift[0],lat+shift[1],u,v,transform=ccrs.PlateCarree(),
						 minlength=10,width=0.005,headlength=5,pivot='middle')
	
		
def drawText(region,ax):
	if region is 'ak':
		ax.text(np.mean(xflux)+8,np.mean(yflux)-1, '   Local\nResource',
						transform=ccrs.Geodetic(), fontsize=12)
		ax.text(np.mean(xflux)+10,np.mean(yflux)-13.5, 'Remote Resource',
						transform=ccrs.Geodetic(), fontsize=12)
	elif region is 'at':
		ax.text(np.mean(xflux)-3.5,np.mean(yflux)+1.5, 'Local Resource',
						transform=ccrs.Geodetic(), fontsize=12)
		ax.text(np.mean(xflux)+18,np.mean(yflux)-2, ' Remote\nResource',
						transform=ccrs.Geodetic(), fontsize=12)
	elif region is 'prusvi':
		ax.text(np.mean(xflux)-1.5,np.mean(yflux)-4.5, 'Local Resource',
						transform=ccrs.Geodetic(), fontsize=12)
		ax.text(np.mean(xflux)+1.5,np.mean(yflux)-7.5, ' Remote\nResource',
						transform=ccrs.Geodetic(), fontsize=12)
	elif region is 'wc':
		ax.text(np.mean(xflux)-2.25,np.mean(yflux)+4.5, '   Local\nResource',
						transform=ccrs.Geodetic(), fontsize=9)
		ax.text(np.mean(xflux)-7,np.mean(yflux)-8.5, 'Remote Resource',
						transform=ccrs.Geodetic(), fontsize=12)
		
			
def padValues(region,x,y,fluxCut=False):
	
	if region is 'prusvi' and fluxCut:
		return x[:-2],y[:-2]
	else:
		if region is 'ak' and not fluxCut:
			x, y = [x.max()+20]+x.tolist(), [y.min()+5]+y.tolist()
			return np.array(x), np.array(y)
		elif region is 'at' and not fluxCut:
			x, y = x.tolist()+[x.min()-30], y.tolist()+[y.min()+10]
			return np.array(x), np.array(y)
		else:
			return x, y
			
			
if __name__ == "__main__":

	
	Text = True
	Legend = False
	Vectors = False
	localArrows = True
	linearAxes = False
	contours = True
	save = False
	
	redAlpha = 0.3
	pc = ccrs.PlateCarree()

	regions = ['ak','at','hi','prusvi','wc']
	regions = ['ak']
		
	for region in regions:
		
		rinf = wr.RegionInfo(region)
		
		extent = rinf.proj.lonlim+rinf.proj.latlim

		fig = plt.figure(10)
		fig.clf()
		
		if linearAxes:
			ax = plt.axes(projection=pc)
		else:
			ax = plt.axes(projection=rinf.proj)
		
		setExtent(extent,region,ax)
		
		inds_all = rinf.con_defs['EEZ']
		for inds in inds_all:
			xflux, yflux = rinf.gridlonlat[:, inds]
			xflux, yflux = padValues(region,xflux,yflux)
			flux = sgeom.LineString(zip(xflux,yflux))
			ax.add_geometries([flux],pc,facecolor='red',
							  	edgecolor='none',alpha=redAlpha)
			
		if region is 'hi': 
			eez = '200'
		else:
			eez = 'eez'
		inds_all = rinf.con_defs[eez]
		for inds in inds_all:
			xflux, yflux = rinf.gridlonlat[:, inds]
			flux = sgeom.LineString(zip(xflux,yflux))
			ax.add_geometries([flux],pc,facecolor='none',
							  	edgecolor='black',linestyle='-',
							    linewidth=2,alpha=1)
			
			
		inds_all = rinf.con_defs['borders']
		
		for inds in inds_all:
			xBound, yBound = rinf.gridlonlat[:, inds]
			xBound, yBound = padValues(region,xBound,yBound,fluxCut=True)
			bound = sgeom.LineString(zip(xBound,yBound))
			ax.add_geometries([bound],pc,facecolor='none',
							  	edgecolor='black',linestyle='--',
							    linewidth=1,alpha=1)
		
			
		sp = cfeature.NaturalEarthFeature(
				category='cultural',
				name='admin_1_states_provinces_lines',
				scale='50m',
				facecolor='none'
			)
		land_50m = cfeature.NaturalEarthFeature(
				'physical','land','50m',
				facecolor=cfeature.COLORS['land'],
				edgecolor='black'
			)
		
		ax.coastlines(resolution='50m',color='black',linewidth=1)
		ax.add_feature(land_50m)
		#ax.add_feature(water_50m)
		ax.add_feature(cfeature.LAKES)
		ax.add_feature(cfeature.OCEAN)
		ax.add_feature(cfeature.RIVERS)
		ax.add_feature(cfeature.BORDERS)
		ax.add_feature(sp,edgecolor='grey')
		
		if contours:
			
			ex = setExtent(extent,region,ax,ret=True)
			
			alpha = 1
			colors = ['black']
			linewidth = 0.5
			linestyle = 'dotted'
			fontsize = 6
			inline = False
			
			levels = {'ak':[[6],[5]],
					 'at':[[6],[5]],
					 'hi':[[6],[5]],
					 'prusvi':[[6],[5]],
					 'wc':[[6],[5]]}
			
			'''
			levels = {'ak':[[37],[32]],
					 'at':[[37],[32]],
					 'hi':[[37],[32]],
					 'prusvi':[[74],[64]],
					 'wc':[[57],[42]]}
			'''
			
			lon,lat = (np.linspace(ex[0],ex[1],200),
					   np.linspace(ex[2],ex[3],100))
			#lon,lat = (np.linspace(-179,181,360),
			#		   np.linspace(-90,90,180))
			
			Lon,Lat = np.meshgrid(lon,lat)
			LonC = ax.contour(Lon,Lat, Lon, levels[region][0][0], transform = pc,
					   			colors=colors,alpha=alpha,
							 	linestyles=linestyle,linewidths=linewidth)
			
			LatC = ax.contour(Lon,Lat, Lat, levels[region][1][0], transform = pc,
							  	colors=colors,alpha=alpha,
							 	linestyles=linestyle,linewidths=linewidth)
			
			ax.clabel(LonC,inline=inline,fontsize=fontsize)
			ax.clabel(LatC,inline=inline,fontsize=fontsize)

		if Vectors:
			vectors(region,extent,ax,localArrows=localArrows)
		
		if Text:
			drawText(region,ax)
		
		if Legend:
		
			local = mpatches.Rectangle((0, 0), 1, 1, facecolor="red",alpha=redAlpha)
			remote = mpatches.Rectangle((0, 0), 1, 1, facecolor=cfeature.COLORS['water'])
			labels = ['Local\nResource',
					  'Remote\nResource']
			plt.legend([local, remote], labels,
					   loc='lower left', ncol=2, bbox_to_anchor=(0.025, -0.1), fancybox=True)

		if linearAxes:
			ax.set_xticks(np.linspace(extent[0], extent[1], 5), crs=pc)
			ax.set_yticks(np.linspace(extent[2], extent[3], 5), crs=pc)
			gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=0.5, color='black', alpha=0.5, linestyle='--')
		
		
		if region is not 'hi' and save is True:
			plt.savefig('./fig/'+region+'_EEZ_resourcePlot.png')	
		
		### Tick labels doesn't work with AlbeirsEqualArea Projection
		#gl.xlabels_top = True
		#gl.ylabels_left = True
		#gl.xlines = False
		#gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
		#gl.xformatter = LONGITUDE_FORMATTER
		#gl.yformatter = LATITUDE_FORMATTER
		#gl.xlabel_style = {'size': 15, 'color': 'gray'}
		#gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
			
		plt.show()
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		