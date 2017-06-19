
package uk.ac.birmingham;

import java.util.Arrays;

import org.scijava.ItemIO;
import org.scijava.command.Command;
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
import net.imglib2.outofbounds.OutOfBoundsMirrorFactory;
import net.imglib2.outofbounds.OutOfBoundsMirrorFactory.Boundary;
import net.imglib2.realtransform.AffineGet;
import net.imglib2.realtransform.AffineRandomAccessible;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.RealViews;
import net.imglib2.type.numeric.RealType;
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
	private OpService opService;

	@Parameter
	private DatasetService datasetService;

	@Parameter(type = ItemIO.OUTPUT)
	private IntervalView<T> bounded;

	@Override
	public void run() {
		final Img<T> image = (Img<T>) currentData.getImgPlus();

		// shear the image
		final AffineTransform3D affine = new AffineTransform3D();
		// final AffineTransform2D affine = new AffineTransform2D();
		affine.set(new double[][] { { 1., 0., 100., 0. }, { 0., 1., 0., 0. }, { 0., 0., 1., 0. }, { 0., 0., 0., 1. } });
		// affine.set(new double[][] { { 1., 1., 0.}, { 0., 1., 0.}, { 0., 0.,
		// 1.} });

		bounded = affineBoundingBox(image, affine, new NLinearInterpolatorFactory<>());
		
		//long[] borderSize = new long[bounded.numDimensions()];
		//Arrays.fill(borderSize, 2);
		
		//RandomAccessibleInterval<FloatType> deconvolved = opService.deconvolve().richardsonLucy(bounded, bounded, borderSize, new OutOfBoundsMirrorFactory<T, RandomAccessibleInterval<T>>(Boundary.DOUBLE), new OutOfBoundsMirrorFactory<T, RandomAccessibleInterval<T>>(Boundary.DOUBLE), new FloatType(), fftType, maxIterations, nonCirculant, accelerate)
		//opService.deconvolve().richardsonLucy(bounded, bounded, 2);
		
	}

	private IntervalView<T> affineBoundingBox(final Img<T> image, final AffineGet affine,
			InterpolatorFactory<T, RandomAccessible<T>> interpolator) {

		
		final RealRandomAccessible<T> field = Views.interpolate(Views.extendZero(image), interpolator);
		final AffineRandomAccessible<T, AffineGet> sheared = RealViews.affine(field, affine);
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
		final Dataset dataset = ij.scifio().datasetIO().open("http://imagej.net/images/FluorescentCells.jpg");

		// show the image
		ij.ui().show(dataset);

		// invoke the plugin
		ij.command().run(LLSPreprocess.class, true);

	}

}
