
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

	@Parameter(label = "Select the image file for the PSF ", persist = true, style = "FILE")
	private File psfFile;

	@Parameter(label = "Define the shear for each slice (pixels): ", persist = false, min = "0.", max = "100.")
	private double shearPix = 1.;

	@Parameter(label = "Define the rotation angle (degrees): ", persist = false, min = "0.", max = "90.")
	private double angleDegrees = 32.8;

	@Parameter(label = "Maximum number of itearions foe deconvolution: ", persist = false, min = "1", max = "200")
	private int maxItDecon = 1;

	@Parameter(type = ItemIO.OUTPUT)
	private IntervalView<FloatType> result;

	@Override
	public void run() {
		final Img<T> image = (Img<T>) currentData.getImgPlus();

		// create affine transform for shear
		final AffineTransform3D affineShear = new AffineTransform3D();
		affineShear.set(new double[][] { { 1., 0., shearPix, 0. }, { 0., 1., 0., 0. }, { 0., 0., 1., 0. },
				{ 0., 0., 0., 1. } });
		
		// apply affine transform and retrieve bounded view
		IntervalView<T> deskewed = affineBoundingBox(image, affineShear, new NLinearInterpolatorFactory<>());
		
		
		try {
			// open the psf stack
			Img<T> psf = (Img<T>) ioService.open(psfFile.getAbsolutePath());
			
			// convert both deskewed View and psf to Float
			// this seems to be needed to get the deconvolution op to work????
			Img<FloatType> deskewedFloat = opService.convert().float32(deskewed);
			Img<FloatType> psfFloat = opService.convert().float32(psf);
			
			// deconvolve with standard RL algorithm
			RandomAccessibleInterval<FloatType> deconvolved = (RandomAccessibleInterval<FloatType>) opService
					.run("richardsonLucy", deskewedFloat, psfFloat, maxItDecon);
			
			// convert rotation angle to radians
			final double angleRad = Math.toRadians(angleDegrees);
		
			// create affine transform for rotation
			final AffineTransform3D affineRot = new AffineTransform3D();
			affineRot.set(new double[][] { { angleRad, 0., 0.5, 0. }, { 0., 1., 0., 0. }, { -0.5, 0., angleRad, 0. },
					{ 0., 0., 0., 1. } });
			// apply affine transform and retrieve bounded view
			IntervalView<FloatType> rotated = affineBoundingBox(deconvolved, affineRot,
					new NLinearInterpolatorFactory<>());
			// Redefine coordinate system so that the image can be displayed despite the possible negative Interval
			// This is necessary because of a bug????
			result = Views.zeroMin(rotated);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	
	/**
	 * The function performs an arbitrary affine transform for an arbitrary number of dimensions.
	 * The resulting view is automatically bounded
	 * 
	 * @param image input image
	 * @param affine defines the transformation
	 * @param interpolator defines the type of interpolation
	 * 
	 * @return the transformed data defined on a View with bounded source
	 */
	private <K extends RealType<K>> IntervalView<K> affineBoundingBox(final RandomAccessibleInterval<K> image,
			final AffineGet affine, InterpolatorFactory<K, RandomAccessible<K>> interpolator) {
		
		// extend with zeros
		final RealRandomAccessible<K> field = Views.interpolate(Views.extendZero(image), interpolator);
		
		// transform the data
		final AffineRandomAccessible<K, AffineGet> sheared = RealViews.affine(field, affine);
		
		int numDims = sheared.numDimensions();
		
		// to hold the bounds of the transformed data
		long[] min = new long[numDims];
		long[] max = new long[numDims];
		Arrays.fill(min, Long.MAX_VALUE);
		Arrays.fill(max, Long.MIN_VALUE);
		
		// to hold the corners of the dataset before and after transformation
		final double[] source = new double[numDims];
		double[] target = new double[numDims];
		
		// dataset dimensions
		final long[] dims = new long[numDims];
		image.dimensions(dims);
		
		// for calculating the corder coordinates
		final long[] normDims = new long[numDims];
		Arrays.fill(normDims, 2);
		
		// to hold unnormalised coordinates, eg [1, 0, 0] or [1, 1, 1]
		final long[] normPos = new long[numDims];
		
		// loop over corners, there are 2^(numDimensions)
		for (int i = 0; i < Math.pow(2, numDims); i++) {
			// calculate unnormalised coordinates
			IntervalIndexer.indexToPosition(i, normDims, normPos);
			// normalise coordinate
			for (int d = 0; d < numDims; d++)
				source[d] = dims[d] * normPos[d];
			// apply affine transform to corner
			affine.apply(source, target);
			
			// if corner is either current extremal value in any dimension update bounds
			for (int d = 0; d < min.length; d++) {
				if (target[d] < min[d])
					min[d] = (long) Math.floor(target[d]);

				if (target[d] > max[d])
					max[d] = (long) Math.ceil(target[d]);

			}
		}
		// create interval from extremal values
		FinalInterval bounds = new FinalInterval(min, max);
		
		// return View with bounded source
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
		
		// crop the dataset so it works on my laptop. This is a bit messy...
		FinalInterval interval = FinalInterval.createMinSize(0, 0, 0, dataset.dimension(0), 400, 30);
		RandomAccessibleInterval<ShortType> croppedAcceesible = (RandomAccessibleInterval<ShortType>) ij.op()
				.run("crop", dataset, interval, true);
		dataset = null;
		// show the cropped image
		ij.ui().show(croppedAcceesible);

		// invoke the plugin
		ij.command().run(LLSPreprocess.class, true);

	}
	


}
