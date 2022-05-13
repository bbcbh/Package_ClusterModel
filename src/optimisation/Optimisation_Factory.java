package optimisation;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.InvalidPropertiesFormatException;
import java.util.Properties;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import sim.SimulationInterface;
import sim.Simulation_ClusterModelGeneration;

public class Optimisation_Factory {

	public static void partnerDistributionFit_Simplex(File baseDir, final double[] init_val_bounded,
			final double[][] boundaries) throws FileNotFoundException, IOException, InvalidPropertiesFormatException {
		
		final File propFile = new File(baseDir, SimulationInterface.FILENAME_PROP);
		final File outputFile = new File(baseDir, String.format("Opt_output_%d.txt", System.currentTimeMillis()));
	
		if (propFile.exists()) {
	
			FileInputStream fIS = new FileInputStream(propFile);
			Properties prop = new Properties();
			prop.loadFromXML(fIS);
			fIS.close();
	
			System.out.println(String.format("Properties file < %s > loaded.", propFile.getAbsolutePath()));
	
			Simulation_ClusterModelGeneration sim = new Simulation_ClusterModelGeneration();
			sim.setBaseDir(baseDir);
			sim.loadProperties(prop);
	
			final long seed;
			final int numSnap;
			final int snapFreq;
			final int[] sampleRange;
	
			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED])) {
				seed = Long
						.parseLong(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));
			} else {
				RandomGenerator rngBase = new MersenneTwisterRandomGenerator(System.currentTimeMillis());
				seed = rngBase.nextLong();
			}
			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])) {
				numSnap = Integer
						.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP]));
			} else {
				numSnap = 1;
			}
			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])) {
				snapFreq = Integer
						.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ]));
			} else {
				snapFreq = AbstractIndividualInterface.ONE_YEAR_INT;
			}
			if (prop.containsKey("POP_PROP_INIT_PREFIX_12")) {
				sampleRange = (int[]) util.PropValUtils.propStrToObject(
						prop.getProperty("POP_PROP_INIT_PREFIX_12"),
						int[].class);
			} else {
				sampleRange = new int[] { 0, snapFreq * numSnap };
			}
	
			MultivariateFunction func = new MultivariateFunction() {
	
				@Override
				public double value(double[] point) {
	
					long tic = System.currentTimeMillis();
	
					Population_Bridging population = new Population_Bridging(seed);
					for (int f = 0; f < Population_Bridging.LENGTH_FIELDS_BRIDGING_POP; f++) {
						if (sim.simFields[f] != null) {
							population.getFields()[f] = sim.simFields[f];
						}
					}
					population.setPrintStatus(null);					
	
					HashMap<String, Object> stepwiseOutput = new HashMap<>();
	
					double[] var_by_gender = new double[Population_Bridging.LENGTH_GENDER];
					Arrays.fill(var_by_gender, 1);
					switch (point.length) {
					case 1:
						// All
						for (int i = 0; i < var_by_gender.length; i++) {
							var_by_gender[i] = point[0];
						}
						break;
					case 2:
						// Hetro female and all msm
						var_by_gender[Population_Bridging.GENDER_HETRO_FEMALE] = point[0];
						var_by_gender[Population_Bridging.GENDER_HETRO_MALE] = point[0];
						var_by_gender[Population_Bridging.GENDER_MSMO] = point[1];
						var_by_gender[Population_Bridging.GENDER_MSMW] = point[1];
						break;
					case 3:
						// Hetro female, MSMO and MSMW
						var_by_gender[Population_Bridging.GENDER_HETRO_FEMALE] = point[0];
						var_by_gender[Population_Bridging.GENDER_HETRO_MALE] = point[0];
						var_by_gender[Population_Bridging.GENDER_MSMO] = point[1];
						var_by_gender[Population_Bridging.GENDER_MSMW] = point[2];
						break;
					default:
						var_by_gender = Arrays.copyOf(point, Population_Bridging.LENGTH_GENDER);
					}
					
					int[] numInPop = (int[]) population.getFields()[Population_Bridging.FIELD_POP_COMPOSITION];
					float[] field_mean_num_partners = (float[]) population
							.getFields()[Population_Bridging.FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS];
					float[] default_mean_num_partners = Arrays.copyOf(field_mean_num_partners, field_mean_num_partners.length);					
					
					int numCat = field_mean_num_partners.length / (1 + Population_Bridging.LENGTH_GENDER);
	
					if (field_mean_num_partners.length == Population_Bridging.LENGTH_GENDER) {
						for (int g = 0; g < field_mean_num_partners.length; g++) {
							field_mean_num_partners[g] = Math.round(field_mean_num_partners[g] * (float) var_by_gender[g]);
						}
					} else {
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							for (int c = 0; c < numCat; c++) {
								int pdIndex = numCat + g * numCat + c;
								field_mean_num_partners[pdIndex] = (float) (field_mean_num_partners[pdIndex] * var_by_gender[g]);
							}
						}
					}
	
					//population.setPrintStatus(System.out);
					population.initialise();
	
					int[] population_num_partner_in_last_12_months_total = new int[field_mean_num_partners.length];
					int sampleCount = 0;
					int startTime = Math.max(0, sampleRange[0] -1);
	
					for (int s = 0; s < numSnap; s++) {
						for (int f = 0; f < snapFreq; f++) {
							if (population.getGlobalTime() == startTime ) {
								population.setStepwise_output(stepwiseOutput);
							} else if (population.getGlobalTime() >= sampleRange[1]) {
								population.setStepwise_output(null);
							}
	
							population.advanceTimeStep(1);
	
							if (population.getGlobalTime() >= sampleRange[0]
									&& population.getGlobalTime() <= sampleRange[1]) {
								int[] population_num_partner_in_last_12_months = (int[]) stepwiseOutput
										.get(Population_Bridging.STEPWISE_OUTPUT_NUM_PARTNERS_IN_12_MONTHS);
								for (int i = numCat; i < population_num_partner_in_last_12_months.length; i++) {
									population_num_partner_in_last_12_months_total[i] += population_num_partner_in_last_12_months[i];
								}
								sampleCount++;
							}
						}
					}
	
					double residue = 0;
					double diff;
					double[] gender_weight = new double[Population_Bridging.LENGTH_GENDER];
					Arrays.fill(gender_weight, 1);
					gender_weight[Population_Bridging.GENDER_HETRO_MALE] = 0;
	
					for (int g = 0; g < gender_weight.length; g++) {
						if (field_mean_num_partners.length == Population_Bridging.LENGTH_GENDER) {
							diff = gender_weight[g] * (((double) population_num_partner_in_last_12_months_total[g]) / sampleCount
									- Math.round(default_mean_num_partners[g] * numInPop[g]));
							residue += Math.pow(diff, 2);
						} else {
							for (int c = 0; c < numCat; c++) {
								int pdIndex = numCat + g * numCat + c;
								diff = gender_weight[g] * (((double)
										population_num_partner_in_last_12_months_total[pdIndex]) / sampleCount
										- Math.round(default_mean_num_partners[pdIndex] * numInPop[g]));
								residue += Math.pow(diff, 2);
							}
						}
					}
					
					String outputString = String.format("x = %s, R = %.4f, Time req. = %.3fs", Arrays.toString(point),
							residue, (System.currentTimeMillis() - tic) / 1000f);
	
					System.out.println(outputString);					
					try {
						PrintWriter pWri = new PrintWriter(new FileWriter(outputFile, true));
						pWri.println(outputString);
						pWri.close();											
					
					}catch(IOException ex) {
						ex.printStackTrace(System.err);
					}
					
					
	
					return residue;
				}
	
			};
	
			MultivariateFunctionMappingAdapter wrapper = new MultivariateFunctionMappingAdapter(func, boundaries[0],
					boundaries[1]);
			SimplexOptimizer optimizer = new SimplexOptimizer(1e-5, 1e-10);
			NelderMeadSimplex simplex = new NelderMeadSimplex(init_val_bounded.length);
			ObjectiveFunction objFunc = new ObjectiveFunction(wrapper);
	
			PointValuePair var = optimizer.optimize(new MaxEval(100), objFunc, simplex, GoalType.MINIMIZE,
					new InitialGuess(wrapper.boundedToUnbounded(init_val_bounded)));
	
			
			String outputString = "Optimised value = " + Arrays.toString(wrapper.unboundedToBounded(var.getPoint()));
			System.out.println(outputString);						
			try {
				PrintWriter pWri = new PrintWriter(new FileWriter(outputFile, true));
				pWri.println(outputString);
				pWri.close();											
			
			}catch(IOException ex) {
				ex.printStackTrace(System.err);
			}
	
		}
	}

}
