
package uk.ac.birmingham;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.io.IOService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccessible;
import net.imglib2.img.Img;
import net.imglib2.interpolation.InterpolatorFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.AffineGet;
import net.imglib2.realtransform.AffineRandomAccessible;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.RealViews;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.ShortType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.IntervalIndexer;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

/**
 * This is a minimal ImageJ command implementing a difference of Gaussian
 * filter.
 */

@Plugin(type = Command.class, menuPath = "Plugins>LLS preprocess")
public class LLSPreprocess<T extends RealType<T>> implements Command {

	@Parameter
	private Dataset currentData;

	@Parameter
	private UIService uiService;

	@Parameter
	private IOService ioService;

	@Parameter
	private OpService opService;

	@Parameter
	private DatasetService datasetService;

	@Parameter(label = "Select the image file for the PSF ", persist = true, style = "FILE")
	private File psfFile;

	@Parameter(type = ItemIO.OUTPUT)
	private IntervalView<FloatType> result;

	@Override
	public void run() {
		final Img<T> image = (Img<T>) currentData.getImgPlus();

		// shear the image
		final AffineTransform3D affine = new AffineTransform3D();
		// final AffineTransform2D affine = new AffineTransform2D();
		affine.set(new double[][] { { 1., 0., 1., 0. }, { 0., 1., 0., 0. }, { 0., 0., 1., 0. }, { 0., 0., 0., 1. } });
		// affine.set(new double[][] { { 1., 1., 0.}, { 0., 1., 0.}, { 0., 0.,
		// 1.} });

		IntervalView<T> bounded = affineBoundingBox(image, affine, new NLinearInterpolatorFactory<>());
		// long[] borderSize = new long[bounded.numDimensions()];
		// Arrays.fill(borderSize, 2);

		try {
			Img<T> psf = (Img<T>) ioService.open(psfFile.getAbsolutePath());
			uiService.show(psf);

			// deconvolved = opService.deconvolve().richardsonLucy(bounded, psf,
			// 5);
			Img<FloatType> boundedFloat = opService.convert().float32(bounded);
			Img<FloatType> psfFloat = opService.convert().float32(psf);
			RandomAccessibleInterval<FloatType> deconvolved = (RandomAccessibleInterval<FloatType>) opService
					.run("richardsonLucy", boundedFloat, psfFloat, 1);
			uiService.show(deconvolved);

			final AffineTransform3D affineRot = new AffineTransform3D();
			affineRot.set(new double[][] { { 0.8660254038, 0., 0.5, 0. }, { 0., 1., 0., 0. },
					{ -0.5, 0., 0.8660254038, 0. }, { 0., 0., 0., 1. } });
			IntervalView<FloatType> rotated = affineBoundingBox(deconvolved, affineRot,
					new NLinearInterpolatorFactory<>());
			result = Views.zeroMin(rotated);
			// uiService.show(rotated);

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private <K extends RealType<K>> IntervalView<K> affineBoundingBox(final RandomAccessibleInterval<K> image,
			final AffineGet affine, InterpolatorFactory<K, RandomAccessible<K>> interpolator) {

		final RealRandomAccessible<K> field = Views.interpolate(Views.extendZero(image), interpolator);
		final AffineRandomAccessible<K, AffineGet> sheared = RealViews.affine(field, affine);
		int numDims = sheared.numDimensions();
		long[] min = new long[numDims];
		long[] max = new long[numDims];
		Arrays.fill(min, Long.MAX_VALUE);
		Arrays.fill(max, Long.MIN_VALUE);

		double[] target = new double[numDims];

		final long[] dims = new long[numDims];
		image.dimensions(dims);

		final long[] normDims = new long[numDims];
		Arrays.fill(normDims, 2);
		final long[] normPos = new long[numDims];
		final double[] source = new double[numDims];
		for (int i = 0; i < Math.pow(2, numDims); i++) {
			IntervalIndexer.indexToPosition(i, normDims, normPos);
			for (int d = 0; d < numDims; d++)
				source[d] = dims[d] * normPos[d];

			affine.apply(source, target);

			for (int d = 0; d < min.length; d++) {
				if (target[d] < min[d])
					min[d] = (long) Math.floor(target[d]);

				if (target[d] > max[d])
					max[d] = (long) Math.ceil(target[d]);

			}
		}
		for (int d = 0; d < min.length; d++) {
			System.out.println(d + " min:" + min[d] + " max:" + max[d]);
		}
		FinalInterval bounds = new FinalInterval(min, max);
		return Views.interval(sheared, bounds);

	}

	/**
	 * This main function serves for development purposes.
	 *
	 * @param args
	 *            whatever, it's ignored
	 * @throws Exception
	 */
	public static void main(final String... args) throws Exception {

		// create the ImageJ application context
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// load example dataset
		Dataset dataset = (Dataset) ij.io().open("testStack.tif");

		FinalInterval interval = FinalInterval.createMinSize(0, 0, 0, dataset.dimension(0), 400, 30);
		RandomAccessibleInterval<ShortType> croppedAcceesible = (RandomAccessibleInterval<ShortType>) ij.op()
				.run("crop", dataset, interval, true);
		dataset = null;
		// show the image
		ij.ui().show(croppedAcceesible);

		// invoke the plugin
		ij.command().run(LLSPreprocess.class, true);

	}

}
