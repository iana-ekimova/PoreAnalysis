import shutil

import numpy
from scipy.io import (loadmat)
#from mat73 import loadmat
#import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import os

import pypore3d
from pypore3d import *
from pypore3d.p3dFiltPy import *
from pypore3d.p3dSkelPy import *
from pypore3d.p3dBlobPy import *
from pypore3d.p3dSITKPy import *
from pypore3d.p3dFiltPy_16 import *
from pypore3d.p3dSITKPy_16 import *

import SimpleITK as sitk
import pypore3d.p3d_SITK_common_lib
import pypore3d.p3d_SITK_read_raw
from pypore3d.p3d_SITK_read_raw import *
from matplotlib.colors import LinearSegmentedColormap



def print_BasicStats(data):
    print("Vv: {}, Sv: {}, Mv: {}, Cv: {}".format(data.Vv, data.Sv, data.Mv, data.Cv))
def print_AnisotropyStats(data):
    print("E: {}, I: {}". format(data.E, data.I))

# https://matplotlib.org/stable/gallery/event_handling/image_slices_viewer.html
class IndexTracker():
    def __init__(self, ax, X, pore_labels=None):
        base_cmap = plt.get_cmap('viridis')
        colors = base_cmap(np.linspace(0, 1, 256 * 256 * 256))
        colors[0] = [0, 0, 0, 1]  # Set the first color to black (R, G, B, A)
        self.custom_cmap = LinearSegmentedColormap.from_list('custom_viridis', colors)

        self.index = 0
        self.X = X
        self.ax = ax
        self.pore_labels = pore_labels
        self.im = ax.imshow(self.X[:, :, self.index], cmap=self.custom_cmap)
        self.vmin = np.min(X)
        self.vmax = np.max(X)
        self.update()

    def on_move(self, event):
        try:
            self.ax.set_title(
                f'Use scroll wheel to navigate\nindex {self.index}, label {self.X[int(event.ydata), int(event.xdata), self.index] - 1000}')
            self.im.axes.figure.canvas.draw()
        except:
            pass

    def on_scroll(self, event):
        increment = 1 if event.button == 'up' else -1
        max_index = self.X.shape[-1] - 1
        self.index = np.clip(self.index + increment, 0, max_index)
        self.update()

    def update(self):
        if self.pore_labels is not None:
            self.ax.clear()
            self.im = self.ax.imshow(self.X[:, :, self.index], cmap=self.custom_cmap, vmin=self.vmin, vmax=self.vmax)
            for label, centroid in self.pore_labels.items():
                z, y, x = centroid
                if label + 1000 == self.X[x, y, self.index]:
                    self.ax.text(y, x, str(label), color='red', fontsize=8, ha='center', va='center')

        self.im.set_data(self.X[:, :, self.index])
        self.im.set_cmap(self.custom_cmap)
        self.ax.set_title(
            f'Use scroll wheel to navigate\nindex {self.index}')
        self.im.axes.figure.canvas.draw()
        self.im.set_clim(self.vmin, self.vmax)


def blob_analysis(small_mat, output):
    pass

FILE_NAME = r"../data/pictures/nestle_5_m_small.mat"
VAR_NAME = ('nestle_5_small')
IMAGE_SPACING = (0.008, 0.008, 0.008)
DIM_X = 400
DIM_Y = 400
DIM_Z = 400

#FILE_NAME = r"../data/pictures/3D_GardenGourmet_21042022_first_small.mat"
#VAR_NAME = "grayscale_image"

pic_small = loadmat(FILE_NAME)
print(pic_small.keys())
small_mat = pic_small[VAR_NAME]

figure_original = plt.figure('Original')
tracker_original = IndexTracker(plt.axes(), small_mat)
figure_original.canvas.mpl_connect('scroll_event', tracker_original.on_scroll)


#plt.show()

"""Filtering small"""
(small_mat * 255).astype("uint8").tofile("tmp.raw")
inImg = py_p3dReadRaw8(r"tmp.raw", small_mat.shape[0], small_mat.shape[1], small_mat.shape[2])
outImg = inImg #outImg = py_p3dMedianFilter8(inImg,400, 400, 400, width = 3 )
#py_p3dWriteRaw8(outImg, r"test.raw", DIM_X, DIM_Y, DIM_Z)
#print("Done")
outImg= py_p3dAutoThresholding8(outImg, DIM_X, DIM_Y, DIM_Z, methodNum = 5)
outImg = py_p3d_Dilate(outImg, DIM_X, DIM_Y, DIM_Z, kWidth = 3)
outImg = py_p3d_Erode(outImg, DIM_X, DIM_Y, DIM_Z, kWidth = 3)

# Store inverted image (will be used later)
invert_vol(outImg, DIM_X, DIM_Y, DIM_Z)
py_p3dWriteRaw8(outImg, r"after_erode.raw", DIM_X, DIM_Y, DIM_Z)
invert_vol(outImg, DIM_X, DIM_Y, DIM_Z)

small_filtered_converted = numpy.fromfile(r"after_erode.raw", dtype=np.uint8) / 256.0
small_filtered_converted = small_filtered_converted.reshape(DIM_X, DIM_Y, DIM_Z)

figure_filtered_small = plt.figure('Filtered 1')
tracker_filtered_small = IndexTracker(plt.axes(), small_filtered_converted)
figure_original.canvas.mpl_connect('scroll_event', tracker_filtered_small.on_scroll)

distImg = outImg
data_output_folder = '/Users/iana/Documents/uni/3d-analysis/scripts'
file_name = "data"
distFile = data_output_folder+file_name+'_distance_field.raw'

"Apply distance field transform"
invert_vol(distImg, DIM_X, DIM_Y, DIM_Z)
outImg16 = py_p3dChamferDT(distImg, DIM_X, DIM_Y, DIM_Z, w1 = 3, w2 = 4, w3 = 5)
invert_vol_16(outImg16, DIM_X, DIM_Y, DIM_Z)
py_p3dWriteRaw16(outImg16, distFile ,DIM_X, DIM_Y, DIM_Z)

distMinFile = data_output_folder+file_name+'_distance_field_min'
py_p3d_HMinimaFilter(distFile, distMinFile, DIM_X, DIM_Y, DIM_Z, threshold = 6)

"Show after HMinimalFilter"
hImg = numpy.fromfile(distMinFile + ".raw", dtype=np.uint16) / 256.0 / 256.0
hImg = hImg.reshape(DIM_X, DIM_Y, DIM_Z)
figure_filtered_hmin = plt.figure('AfterHMinimal')
tracker_filtered_hmin = IndexTracker(plt.axes(), hImg)
figure_original.canvas.mpl_connect('scroll_event', tracker_filtered_hmin.on_scroll)

"Apply Watershed"
watershed_image_path = data_output_folder+file_name+'_watershed'
py_p3d_WatershedSegmentation(distMinFile+".raw", watershed_image_path, DIM_X, DIM_Y, DIM_Z, level = 0, connected = False)

"Show watershed image (stored as 32 bit by SITK)"
wshImg = numpy.fromfile(watershed_image_path + ".raw", dtype=np.uint32)
wshImg = wshImg.reshape(DIM_X, DIM_Y, DIM_Z)
figure_filtered_wsh = plt.figure('Watershed')
tracker_filtered_wsh = IndexTracker(plt.axes(), (wshImg != 0) * 1000 + wshImg)
figure_original.canvas.mpl_connect('scroll_event', tracker_filtered_wsh.on_scroll)


"Apply mask to separate watershed regions"
#invert_vol(image_after_erode, DIM_X, DIM_Y, DIM_Z)
bd_img = Read_Raw("after_erode.raw", [DIM_X, DIM_Y, DIM_Z], image_spacing=IMAGE_SPACING)
wsd_img = Read_Raw(watershed_image_path+'.raw' , [DIM_X, DIM_Y, DIM_Z], sitk.sitkUInt32, image_spacing=IMAGE_SPACING)

wsd_labeled_filter = sitk.RelabelComponentImageFilter()
wsd_labeled_img = wsd_labeled_filter.Execute(wsd_img)

md_img = sitk.Mask(wsd_labeled_img, bd_img)

# Border filter
bgp = sitk.BinaryGrindPeak(md_img != 0)
#md_img = sitk.MaskNegated(md_img, bgp)

md_image_path = data_output_folder+file_name
sitk.WriteImage(md_img, md_image_path+'.mhd')

########

"Show watershed masked image (stored as 32 bit by SITK)"
wshImg = numpy.fromfile(md_image_path + ".raw", dtype=np.uint32)
shutil.copy(md_image_path + ".raw", os.path.join("/Users/iana/Documents/uni/3d-analysis/results", f"{VAR_NAME}_filtered.raw"))
#wshImg = np.array((wshImg == 0) * 255, dtype=np.uint8) # Convert to binary with inversion
wshImg = wshImg.reshape(DIM_X, DIM_Y, DIM_Z)
figure_filtered_wsh2 = plt.figure('Masked')
tracker_filtered_wsh2 = IndexTracker(plt.axes(), (wshImg != 0) * 1000 + wshImg)
figure_original.canvas.mpl_connect('scroll_event', tracker_filtered_wsh2.on_scroll)

#plt.show()

shape_stats = sitk.LabelShapeStatisticsImageFilter()
shape_stats.ComputeOrientedBoundingBoxOn()
shape_stats.ComputeFeretDiameterOn()
shape_stats.Execute(md_img)
#print(shape_stats)

pore_labels = {
    i: md_img.TransformPhysicalPointToIndex(shape_stats.GetCentroid(i))
    for i in shape_stats.GetLabels()
}

print(pore_labels)

#pore_labels = {label: (x, y, z) for label, (x, y, z) in pore_labels.items()}

figure_filtered_wsh2_masked = plt.figure('Masked with Labels')
tracker_filtered_wsh2_masked = IndexTracker(plt.axes(), (wshImg != 0) * 1000 + wshImg, pore_labels=pore_labels)
figure_original.canvas.mpl_connect('scroll_event', tracker_filtered_wsh2_masked.on_scroll)
figure_filtered_wsh2_masked.canvas.mpl_connect('scroll_event', tracker_filtered_wsh2_masked.on_scroll)
figure_filtered_wsh2_masked.canvas.mpl_connect('motion_notify_event', tracker_filtered_wsh2_masked.on_move)

plt.show()

# Define column names
column_names = [
    'PoreNumber',
    'PhysicalSize',
    'Perimeter',
    'Elongation',
    'Flatness',
    'Roundness',
    'Equivalent Ellipsoid Diameter',
    'Equivalent Spherical Radius',
    'Equivalent Spherical Perimeter',
    'OrientedBoundingBoxSize_X',
    'OrientedBoundingBoxSize_Y',
    'OrientedBoundingBoxSize_Z',
    'Centroid_X',
    'Centroid_Y',
    'Centroid_Z'
]

# Prepare the data
stats_list = [
    (
        i,
        shape_stats.GetPhysicalSize(i),
        shape_stats.GetPerimeter(i),
        shape_stats.GetElongation(i),
        shape_stats.GetFlatness(i),
        shape_stats.GetRoundness(i),
        shape_stats.GetEquivalentEllipsoidDiameter(i),
        shape_stats.GetEquivalentSphericalRadius(i),
        shape_stats.GetEquivalentSphericalPerimeter(i),
        shape_stats.GetOrientedBoundingBoxSize(i)[0],
        shape_stats.GetOrientedBoundingBoxSize(i)[1],
        shape_stats.GetOrientedBoundingBoxSize(i)[2],
        md_img.TransformPhysicalPointToIndex(shape_stats.GetCentroid(i))[2],
        md_img.TransformPhysicalPointToIndex(shape_stats.GetCentroid(i))[1],
        md_img.TransformPhysicalPointToIndex(shape_stats.GetCentroid(i))[0],
    )
    for i in shape_stats.GetLabels()
]

df = pd.DataFrame(stats_list, columns=column_names)
df_sorted = df.sort_values(by='PhysicalSize', ascending=True)
df_sorted['CumulativeSum'] = df_sorted['PhysicalSize'].cumsum()
total_physical_size = df_sorted['PhysicalSize'].sum()
df_sorted['CumulativePercentage'] = (df_sorted['CumulativeSum'] / total_physical_size) * 100
df_filtered = df_sorted[df_sorted['CumulativePercentage'] >= 1.5]
#df_filtered = df_filtered.drop(columns=['CumulativeSum', 'CumulativePercentage'])

print(df_filtered)

filename = f"{VAR_NAME}_fss.csv"
output_path = os.path.join("/Users/iana/Documents/uni/3d-analysis/results", filename)
df_filtered.to_csv(output_path, index=False)

print(f"Data has been written to {output_path}")
plt.show()
exit()
