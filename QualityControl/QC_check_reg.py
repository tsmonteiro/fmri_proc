#!/usr/bin/env python



import argparse


import nibabel as nib


import matplotlib.pyplot as plt
import numpy as np


import nilearn.plotting as nlp
import nilearn.image as nimg


from PIL import Image, ImageDraw

def save_normalisation_gif(imgFile1, imgFile2, outFile, thr_top=200, thr_base=0):
	#baseDir = 'C:/Users/u0101486/Documents/workspace/tmp/CRSN_084'
	imBaseFile  =  imgFile1
	imTopFile   =  imgFile2


	imBase = nib.load(imBaseFile)
	imBase = imBase.get_data()
	if len(imBase.shape) > 3:
		imBase = imBase[:,:,:,0]


	imTop = nib.load(imTopFile)
	imTop = imTop.get_data()

	if len(imTop.shape) > 3:
		imTop = imTop[:,:,:,0]

	# Use this to threshold image
	np.where(imTop<thr_top, 0, imTop)
	np.where(imBase<thr_base, 0, imBase)


	x,y,z = imBase.shape

	mZ = int(np.round(z/2))
	mY = int(np.round(y/2))
	mX = int(np.round(x/2))

	imBase1 = np.concatenate( (imBase[:,mY-20,:],
		                  imBase[:,mY-10,:],
		                  imBase[:,mY,:],
		                  imBase[:,mY+10,:],
		                  imBase[:,mY+20,:],), axis=1 )

	imBase2 = np.concatenate( (imBase[mX-20,:,:],
		                  imBase[mX-10,:,:],
		                  imBase[mX-0,:,:],
		                  imBase[mX+10,:,:],
		                  imBase[mX+20,:,:],), axis=1 )

	imBase = np.concatenate( (imBase1, imBase2 ), axis=0 )
	del imBase1
	del imBase2


	imTop1 = np.concatenate( (imTop[:,mY-20,:],
		                  imTop[:,mY-10,:],
		                  imTop[:,mY,:],
		                  imTop[:,mY+10,:],
		                  imTop[:,mY+20,:],), axis=1 )

	imTop2 = np.concatenate( (imTop[mX-20,:,:],
		                  imTop[mX-10,:,:],
		                  imTop[mX-0,:,:],
		                  imTop[mX+10,:,:],
		                  imTop[mX+20,:,:],), axis=1 )

	imTop = np.concatenate( (imTop1, imTop2 ), axis=0 )

	del imTop1
	del imTop2

	maxBase = np.max(np.max(imBase))
	imBase = imBase/maxBase
	imBase = (imBase)*255.


	maxTop = np.max(np.max(imTop))
	imTop = imTop/maxTop
	imTop = (imTop)*255.


	h,w = np.transpose(imTop).shape
	im1 = Image.fromarray( imBase.astype('uint8') ).rotate(90,expand=1)
	im2 = Image.fromarray( imTop.astype('uint8') ).rotate(90,expand=1)



	draw = ImageDraw.Draw(im1)    
	draw.ellipse((0, 0, 20, 20), fill=(255))
	del draw

	draw = ImageDraw.Draw(im2)    
	draw.ellipse((0, 0, 20, 20), fill=(50))
	del draw



	#im = merge_images(imBase, imTop, 0.5)
	alphaInc = np.linspace(0, 1, 10)
	alphaDec = np.linspace(1, 0, 10)
	ims = []
	for alpha in alphaInc:
		ims.append(Image.blend(im1,im2,alpha).resize((w*2,h*2), resample=Image.BILINEAR)   )


	ims.append(Image.blend(im1,im2,1).resize((w*2,h*2), resample=Image.BILINEAR)  )
	ims.append(Image.blend(im1,im2,1).resize((w*2,h*2), resample=Image.BILINEAR)  )
	ims.append(Image.blend(im1,im2,1).resize((w*2,h*2), resample=Image.BILINEAR)  )
	ims.append(Image.blend(im1,im2,1).resize((w*2,h*2), resample=Image.BILINEAR)  )
	ims.append(Image.blend(im1,im2,1).resize((w*2,h*2), resample=Image.BILINEAR)  )

	for alpha in alphaDec:
	    ims.append(Image.blend(im1,im2,alpha).resize((w*2,h*2), resample=Image.BILINEAR)  )


	ims.append(Image.blend(im1,im2,0).resize((w*2,h*2), resample=Image.BILINEAR)  )
	ims.append(Image.blend(im1,im2,0).resize((w*2,h*2), resample=Image.BILINEAR)  )
	ims.append(Image.blend(im1,im2,0).resize((w*2,h*2), resample=Image.BILINEAR)  )
	ims.append(Image.blend(im1,im2,0).resize((w*2,h*2), resample=Image.BILINEAR)  )
	ims.append(Image.blend(im1,im2,0).resize((w*2,h*2), resample=Image.BILINEAR)  )

	ims[0].save(outFile, format='GIF', 
	   append_images=ims[1:], save_all=True, duration=100, loop=0)


def save_normalisation_plot(imgPath, edgePath, outFile, figDpi):

	baseImg = nib.load(imgPath)
	data    = np.array(baseImg.get_data())



	if len(data.shape) > 3:
		data = np.mean( data, 3)

	vmax = np.amax(data)

	vmin = vmax * 0
	vmax = vmax * 0.9


	zCuts=(-34,58)
	xCuts=(-40,40)
	yCuts=(-80,40)

	#TODO [10.9.19] If not in MNI space, cuts likely differ and might need to be taken into account
	
	niimg = nimg.new_img_like(baseImg, data)


	fig = plt.figure(figsize=(30,6), dpi=figDpi, facecolor='w', edgecolor='k')

	# Display the different planes
	cuts = np.linspace(zCuts[0], zCuts[1], 20)

	display1 = nlp.plot_epi(niimg, cmap='gray', figure=fig, cut_coords=cuts, display_mode='z', draw_cross=False, vmax=vmax, vmin=vmin, \
			axes=(0,0,1,.33))
	display1.add_edges(edgePath)

	cuts = np.linspace(xCuts[0], xCuts[1], 20)

	display2 = nlp.plot_epi(niimg, cmap='gray', figure=fig, cut_coords=cuts, display_mode='x', draw_cross=False, vmax=vmax, vmin=vmin, \
			axes=(0,.33,1,.33))
	display2.add_edges(edgePath)


	cuts = np.linspace(yCuts[0], yCuts[1], 20)

	display3 = nlp.plot_epi(niimg, cmap='gray', figure=fig, cut_coords=cuts, display_mode='y', draw_cross=False, vmax=vmax, vmin=vmin, \
			axes=(0,.66,1,.33))
	display3.add_edges(edgePath)

	plt.savefig(outFile)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF function definitions
#
# ++++++++++++++++++++++++++++++++++++++++++++++


parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options                    
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-a', '-in', dest="inDir", required=True, help='Dir where EPI + masks are stored [MNI SPACE]' )

reqoptions.add_argument('-x', '-outf', dest="outFname", required=False, default='greyplot.png', help='Name of the output plot' )


reqoptions.add_argument('-k', '-im1', dest="im1", required=False,  help='Image 1 in normalisation check' )
reqoptions.add_argument('-l', '-im2', dest="im2", required=False,  help='Image 2 in normalisation check' )



reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )


args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir



# PNG resolution of the saved file
figDpi=int(args.dpi)

outFile = outDir + '/' + args.outFname

outGFile = outDir + '/' + args.outFname + '.gif'



# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )




# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF PARAMETER/SETUP
#
# ++++++++++++++++++++++++++++++++++++++++++++++++



save_normalisation_plot(args.im1, args.im2, outFile, figDpi=figDpi)
#save_normalisation_gif(args.im1, args.im2, outGFile, thr_top=10, thr_base=0)


	
