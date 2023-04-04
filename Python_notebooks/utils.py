import pandas as pd
import untangle
import numpy as np
from PIL import Image
from scipy import interpolate
import PIL
from tqdm import tqdm as tqdm

import matplotlib.pyplot as plt
import matplotlib.image

PIL.Image.MAX_IMAGE_PIXELS = 1063733067


def poly_area(x, y):
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


class DVPXML:
    def __init__(self, path):
        self.path = path
        self.content = untangle.parse(path)
        self.parse_shapes()
        self.read_calibration_points()

    def parse_shapes(self):
        self.n_shapes = int(self.content.ImageData.ShapeCount.cdata)

    def return_shape(self, index):
        if index > self.n_shapes:
            raise ValueError(f"Maximum shape is {self.n_shapes}")

        n_points = int(
            self.content.ImageData.__getattr__(f"Shape_{index}").PointCount.cdata
        )
        pts = np.zeros((n_points, 2))

        for i in range(n_points):
            pts[i, 0] = float(
                self.content.ImageData.__getattr__(f"Shape_{index}")
                .__getattr__(f"X_{i+1}")
                .cdata
            )
            pts[i, 1] = float(
                self.content.ImageData.__getattr__(f"Shape_{index}")
                .__getattr__(f"Y_{i+1}")
                .cdata
            )
        return pts[:, 0], pts[:, 1]

    def read_calibration_points(self):
        self.y_calibration = [
            int(self.content.ImageData.__getattr__(f"Y_CalibrationPoint_{i+1}").cdata)
            for i in range(3)
        ]
        self.x_calibration = [
            int(self.content.ImageData.__getattr__(f"X_CalibrationPoint_{i+1}").cdata)
            for i in range(3)
        ]


class DVPMETA:
    def __init__(self, path):
        self.path = path
        self.metadata = pd.read_csv(path, sep="\t")

    def slice_subset(self, selected_slide):
        metadata = self.metadata.copy()
        sub = metadata[metadata["Slide"] == selected_slide]
        return sub


class ImXML:
    def __init__(self, METADATA_PATH, XML_PATH, IM_PATH):

        self.dvpmeta = DVPMETA(METADATA_PATH)
        self.dvpxml = DVPXML(XML_PATH)
        self.im_path = IM_PATH
        self.load_image()

    def load_image(self):

        self.im = np.array(Image.open(self.im_path))
        self.im_shape = self.im.shape

    def bounding_rect(self, x, y):

        return [
            int(np.floor(min(x))),
            int(np.ceil(max(x))),
            int(np.floor(min(y))),
            int(np.ceil(max(y))),
        ]

    def calibration(self, slide):
        self.slide = slide

        sub = self.dvpmeta.slice_subset(slide)
        resolution = sub["resolution"].iloc[0]

        xx = sub["X"] * 1 / resolution + self.im_shape[1] / 2  # this is in pixels
        yy = sub["Y"] * 1 / resolution + self.im_shape[0] / 2

        self.calib_x = xx
        self.calib_y = yy

        self.fxx = interpolate.interp1d(
            self.dvpxml.y_calibration, xx, fill_value="extrapolate"
        )
        self.fyy = interpolate.interp1d(
            self.dvpxml.x_calibration, yy, fill_value="extrapolate"
        )

    def show_all_shapes(self):

        plt.figure(figsize=(10, 10))

        plt.imshow(self.im, vmin=0, vmax=1000, cmap="gray")

        for i in range(1, self.dvpxml.n_shapes):

            try:
                x, y = self.dvpxml.return_shape(i)
            except AttributeError:
                break

            x_px = self.fxx(y)
            y_px = self.fyy(x)

            plt.plot(x_px, y_px, color="blue")

        plt.scatter(self.calib_x, self.calib_y, color="r", marker="+")
        plt.scatter(
            self.fxx(self.dvpxml.y_calibration),
            self.fyy(self.dvpxml.x_calibration),
            color="r",
            marker="x",
        )
        plt.title(self.im_path)
        plt.show()

    def slice_shape(self, idx, mode="bounding", plot=False):

        x, y = self.dvpxml.return_shape(idx)
        x_px = self.fxx(y)
        y_px = self.fyy(x)

        br = self.bounding_rect(x_px, y_px)

        selection = self.im[br[2] : br[3], br[0] : br[1]]

        if plot:
            plt.figure()
            plt.imshow(selection)  # , vmin = 0, vmax = 500)
            plt.plot(x_px - br[0] - 0.5, y_px - br[2] - 0.5, "r", linestyle=":")
            plt.title(f"Shape {idx}")
            plt.show()

        return selection.astype("uint64")

    def slice_shapes(self, mode="bounding"):

        results = []

        for i in tqdm(range(1, self.dvpxml.n_shapes)):
            try:
                selection = self.slice_shape(i)
            except AttributeError:
                break

            x, y = self.dvpxml.return_shape(i)
            area = poly_area(x, y)

            results.append(
                (
                    i,
                    np.mean(selection),
                    np.median(selection),
                    np.min(selection),
                    np.max(selection),
                    area,
                    self.im_path,
                    self.slide,
                )
            )

        return pd.DataFrame(
            results,
            columns=["Index", "Mean", "Median", "Min", "Max", "Area", "Image", "Slide"],
        )

    def slice_shapes_(self, mode="bounding"):

        start = time()

        results = []

        for i in tqdm(range(1, self.dvpxml.n_shapes)):
            try:
                selection = self.slice_shape(i)
            except AttributeError:
                break

            x, y = self.dvpxml.return_shape(i)
            area = poly_area(x, y)

            results.append(
                (
                    i,
                    np.mean(selection),
                    np.median(selection),
                    np.min(selection),
                    np.max(selection),
                    area,
                    self.im_path,
                    self.slide,
                )
            )

        return pd.DataFrame(
            results,
            columns=["Index", "Mean", "Median", "Min", "Max", "Area", "Image", "Slide"],
        )

    def export_shapes(self, path, mode="bounding", offset=0):
        for i in tqdm(range(1, self.dvpxml.n_shapes)):
            try:

                x, y = self.dvpxml.return_shape(i)
                x_px = self.fxx(y)
                y_px = self.fyy(x)

                br = self.bounding_rect(x_px, y_px)

                xmin, xmax = br[2], br[3]
                ymin, ymax = br[0], br[1]

                if offset != 0:

                    x_mean = int((xmin + xmax) / 2)
                    y_mean = int((ymin + ymax) / 2)

                    xmin = x_mean - offset
                    ymin = y_mean - offset

                    xmax = x_mean + offset
                    ymax = y_mean + offset

                if xmin < 0:
                    xmin = 0
                if ymin < 0:
                    ymin = 0

                if xmax > self.im.shape[0]:
                    xmax = self.im.shape[0]

                if ymax > self.im.shape[1]:
                    ymax = self.im.shape[1]

                if offset != 0:

                    cp = self.im.copy()

                    cp[br[2] : br[3], br[0] : br[1]] = self.im[
                        xmin:xmax, ymin:ymax
                    ].max()
                    cp[br[2] + 1 : br[3] - 1, br[0] + 1 : br[1] - 1] = self.im[
                        br[2] + 1 : br[3] - 1, br[0] + 1 : br[1] - 1
                    ]

                    selection = cp[xmin:xmax, ymin:ymax]

                    xmin_ = br[2] - x_mean
                    xmax_ = br[3] - x_mean

                    ymin_ = br[0] - y_mean
                    ymax_ = br[1] - y_mean

                    matplotlib.image.imsave(f"{path}_{i}.png", selection, cmap="gray")
                else:
                    selection = self.im[xmin:xmax, ymin:ymax]

                    plt.figure()
                    plt.imshow(selection)  # , vmin = 0, vmax = 500)
                    plt.plot(x_px - br[0] - 0.5, y_px - br[2] - 0.5, "r", linestyle=":")
                    plt.axis("off")
                    plt.savefig(f"{path}_{i}.png", bbox_inches="tight")
                    # plt.show()
                    plt.close()
            except AttributeError as e:
                print(e)
                break
