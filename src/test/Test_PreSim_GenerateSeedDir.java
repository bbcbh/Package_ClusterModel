package test;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import random.MersenneTwisterRandomGenerator;
import sim.Simulation_ClusterModelGeneration;
import util.Util_7Z_CSV_Entry_Extract_Callable;

public class Test_PreSim_GenerateSeedDir {

	public static void main(String[] args) throws IOException {
		final String REGEX_CMAP = Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP.replaceAll("%d",
				"(-{0,1}(?!0)\\\\d+)");

		File cMapDir, templateDir, sourcePBS, lookup_file = null;
		String dirFormat;
		String seedHeader, paramStr;
		String outputDirName;
		String idLine_str;
		double[][] param_range;
		int numSimPerDir = 16;
		int numDirPerMap = 10;
		int numMapInclMax = Integer.MAX_VALUE;

		final int PROP_TYPE_SIM = 0;
		final int PROP_TYPE_SIM_MSM_DOXY_PEP = PROP_TYPE_SIM + 1;
		final int PROP_TYPE_SIM_MSM_VIABILITY = PROP_TYPE_SIM_MSM_DOXY_PEP + 1;
		final int PROP_TYPE_SIM_BALI = PROP_TYPE_SIM_MSM_VIABILITY + 1;
		final int PROP_TYPE_OPT = PROP_TYPE_SIM_BALI + 1;
		final int PROP_TYPE_OPT_BALI = PROP_TYPE_OPT + 1;
		final int PROP_TYPE_OPT_MSM = PROP_TYPE_OPT_BALI + 1;
		final int PROP_TYPE_OPT_MSM_VIABILITY = PROP_TYPE_OPT_MSM + 1;

		int propType = PROP_TYPE_SIM_MSM_VIABILITY;

		String[] cMap_seedCopy_name = new String[] {}; // new String[] { "RiskGrp_Map_%s.csv" };
		File[] direct_copy = null;

		long rng_seed = 22519122707291119l;
		int extra_seed = 0;

		switch (propType) {
		case PROP_TYPE_OPT_MSM:
			outputDirName = "Opt_Shift_MSM";
			idLine_str = "#PBS -N OptMSM%03d";
			sourcePBS = new File(
					"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Opt_Template.pbs");
			cMapDir = new File("C:\\Users\\Bhui\\OneDrive - UNSW\\Bridging_model\\Test_cMap_MSM");
			templateDir = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Opt_MSM_Multi");
			dirFormat = "Opt_Extra_%03d";
			seedHeader = "CMAP_SEED,SIM_SEED";
			paramStr = "";
			param_range = null;
			numSimPerDir = 2;
			numMapInclMax = 100;
			numDirPerMap = 8;

			break;
		case PROP_TYPE_OPT_MSM_VIABILITY:
			outputDirName = "Opt_Shift_MSM_VIABILITY";
			idLine_str = "#PBS -N OptVAL%03d";
			sourcePBS = new File(
					"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Opt_Template.pbs");
			cMapDir = new File("C:\\Users\\Bhui\\OneDrive - UNSW\\Bridging_model\\Test_cMap_MSM");
			templateDir = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Opt_MSM_Multi_Viability");
			dirFormat = "Opt_Extra_%03d";
			seedHeader = "CMAP_SEED,SIM_SEED,"
					+ "PROP_PROB_NON_VIABLE_TREATMENT_1_1,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_1_1,"
					+ "PROP_PROB_NON_VIABLE_TREATMENT_1_2,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_1_2,"
					+ "PROP_PROB_NON_VIABLE_TREATMENT_1_3,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_1_3,"
					+ "PROP_PROB_NON_VIABLE_TREATMENT_2_1,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_2_1,"
					+ "PROP_PROB_NON_VIABLE_TREATMENT_2_2,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_2_2,"
					+ "PROP_PROB_NON_VIABLE_TREATMENT_2_3,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_2_3,"
					+ "PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION_1_2_0,PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION_1_3_0,"
					+ "PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION_2_2_0,PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION_2_3_0";
			paramStr = ",0.0,1.00,0.0,1.00,0.0,1.00,0.0,1.00,0.0,1.00,0.0,1.00,10,10,10,10";
			param_range = null;
			numSimPerDir = 2;
			numMapInclMax = 100;
			numDirPerMap = 8;
			break;
		case PROP_TYPE_OPT_BALI:
			outputDirName = "Opt_Shift_Bali";
			idLine_str = "#PBS -N OptBA_%03d";
			sourcePBS = new File(
					"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Opt_Template.pbs");
			cMapDir = new File("C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\GenMap_BALI");
			templateDir = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Opt_Bali");
			dirFormat = "Opt_Extra_%03d";
			seedHeader = "CMAP_SEED,SIM_SEED";
			paramStr = "";
			param_range = null;
			numSimPerDir = 2;
			numMapInclMax = 100;
			numDirPerMap = 5;
			cMap_seedCopy_name = new String[0];
			break;
		case PROP_TYPE_OPT:
			outputDirName = "Opt_Shift_ALL";
			idLine_str = "#PBS -N Opt_%03d";
			sourcePBS = new File(
					"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Opt_Template.pbs");
			cMapDir = new File("C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\Test_cMap");
			templateDir = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Opt");
			dirFormat = "Opt_Extra_%03d";
			seedHeader = "CMAP_SEED,SIM_SEED";
			paramStr = "";
			param_range = null;
			numSimPerDir = 2;
			numDirPerMap = 50;
			break;
		case PROP_TYPE_SIM_MSM_DOXY_PEP:
			outputDirName = "MultiTransmission_MSM_S3_PA8_180_075";
			idLine_str = "#PBS -N MSM_S2U2_%02d";
			templateDir = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Sim_DoxyPEP");
			sourcePBS = new File(
					"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Sim_Template.pbs");
			cMapDir = new File("C:\\Users\\Bhui\\OneDrive - UNSW\\Bridging_model\\Test_cMap_MSM");
			lookup_file = null;
			direct_copy = new File[] { new File(templateDir, "simSpecificSim.prop"),
					new File(templateDir, "simSpecificSwitch.prop") };
			// direct_copy = new File[] { new File(templateDir, "simSpecificSim.prop")};

			dirFormat = "SimClusterModel_Transmission_Extra_%03d";
			seedHeader = "CMAP_SEED,SIM_SEED," + "16_1_16,18_4_32,19_31_8," // TP
					+ "16_4_16,16_8_16,16_16_16,16_32_16,16_64_16,16_128_16,16_256_16,18_2048_32,18_16384_32,18_131072_32,19_32_8,19_64_8," // NG
					+ "16_512_16,16_1024_16,16_2048_16,18_1048576_32,18_8388608_32,18_67108864_32"; // CT

			paramStr = "," + "0.019000,450.000000,0.575000" // TP
					+ ",0.549154,0.339017,0.231058,0.066618,0.029145,0.012548,0.003501,61.1412550,154.422144,108.275146,0.796314,0.443171" // NG
					+ ",0.128574,0.049310,0.010380,384.954040,618.176299,289.239346" // CT
					+ "";

			param_range = new double[][] { new double[] { 0.025, -1 } };
			numSimPerDir = 16;
			numDirPerMap = 5;
			break;
		case PROP_TYPE_SIM_MSM_VIABILITY:
			extra_seed = 19;
			int seedOffset = 0;
			outputDirName = String.format("MSM_Viability_Baseline_All_%03d",extra_seed-seedOffset);
			idLine_str = "#PBS -N MSM_VTB"+ (extra_seed-seedOffset)+"_%02d";
			templateDir = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Sim_MSM_Viablity");
			sourcePBS = new File(
					"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Sim_Template.pbs");
			cMapDir = new File("C:\\Users\\Bhui\\OneDrive - UNSW\\Bridging_model\\Test_cMap_MSM");
			lookup_file = null;
			direct_copy = new File[] { new File(templateDir, "simSpecificSim.prop") };
			dirFormat = "MSM_Viability_Baseline_All_000_%03d";

			seedHeader = "CMAP_SEED,SIM_SEED," + "16_1_16,18_4_32,19_31_8" // TP
					+ ",16_4_16,16_8_16,16_16_16,16_32_16,16_64_16,16_128_16,16_256_16" // NG trans
					+ ",18_2048_32,18_16384_32,18_131072_32,19_32_8,19_64_8" // NG dur, sym
					+ ",16_512_16,16_1024_16,16_2048_16" // CT trans
					+ ",18_1048576_32,18_8388608_32,18_67108864_32" // CT dur
					+ ",PROP_DUR_ADJ_NON_VIABLE_TREATMENT_1_2" + ",PROP_DUR_ADJ_NON_VIABLE_TREATMENT_1_3" // NG Non-viability duration
					+ ",PROP_DUR_ADJ_NON_VIABLE_TREATMENT_2_1" + ",PROP_DUR_ADJ_NON_VIABLE_TREATMENT_2_2" + ",PROP_DUR_ADJ_NON_VIABLE_TREATMENT_2_3" // CT Non-viability duration
//					+ ",PROP_PROB_NON_VIABLE_TREATMENT_1_1,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_1_1" // Non-viability from treatment: NG_Genital
//					+ ",PROP_PROB_NON_VIABLE_TREATMENT_1_2,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_1_2" // Non-viability from treatment: NG_Rectal
//					+ ",PROP_PROB_NON_VIABLE_TREATMENT_1_3,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_1_3" // Non-viability from treatment: NG_Orophraygneal
//					+ ",PROP_PROB_NON_VIABLE_TREATMENT_2_1,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_2_1" // Non-viability from treatment: CT_Genital
//					+ ",PROP_PROB_NON_VIABLE_TREATMENT_2_2,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_2_2" // Non-viability from treatment: CT_Rectal
//					+ ",PROP_PROB_NON_VIABLE_TREATMENT_2_3,PROP_DUR_ADJ_NON_VIABLE_TREATMENT_2_3" // Non-viability from treatment: CT_Orophraygneal
//					+ ",PROP_PROB_NON_VIABLE_TRANSMISSION_1_1_2,PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION_1_2_0" // Non-viability from transmission: NG genital -> rectal
//					+ ",PROP_PROB_NON_VIABLE_TRANSMISSION_1_1_3,PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION_1_3_0" // Non-viability from transmission: NG genital -> orophraygneal
//					+ ",PROP_PROB_NON_VIABLE_TRANSMISSION_2_1_2,PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION_2_2_0" // Non-viability from transmission: CT genital -> rectal
//					+ ",PROP_PROB_NON_VIABLE_TRANSMISSION_2_1_3,PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION_2_3_0"  // Non-viability from transmission: CT genital -> orophraygneal		
					+ "";

			// Pre Parameter Sel
//			paramStr = ",0.019,450.0,0.5750" // TP
//					+ ",0.549154,0.339017,0.231058,0.066618,0.029145,0.012548,0.003501" // NG trans
//					+ ",61.1412550,154.422144,108.275146,0.796314,0.443171" // NG dur, sym
//					+ ",0.128574,0.049310,0.010380,384.954040,618.176299,289.239346" // CT
//					+ ",5" // Non-viability from treatment: NG_Rectal
//					+ ",60" // Non-viability from treatment: NG_Orophraygneal
//					+ "";

			paramStr = "0.018486482800664228,467.5382190942764,0.5534875129152894" // TP
					+ ",0.8880828619003296,0.02171908180325488,0.6604948802639541,0.124252223101718" // NG Trans
					+ ",0.3124454112949451,0.06683035078550102,0.003449702262878418" // NG Trans Cont.
					+ ",39.43349003791809,129.6459424495697,147.29536056518555,0.8111539855256751,0.4523859110241748" // NG
																														// dur,
																														// sym
					+ ",0.12858378217022892,0.047125692389754906,0.010213328315341386" // CT trans
					+ ",369.0054189243492,630.1818077835196,280.3020495646033" // CT dur
					+ ",5.274476647377014" // Non-viability from treatment: NG_Rectal
					+ ",49.20526874065399" // Non-viability from treatment: NG_Orophraygneal
					+ ",0.5" // Non-viability from treatment: CT_Urethral
					+ ",0.5" // Non-viability from treatment: CT_Rectal
					+ ",0.5" // Non-viability from treatment: CT_Orophraygneal
					+ "";

			// Usage: Pre-define range:
			// new double[][] { new double[] { min_range, max_range} }
			// or new double[][] { new double[] { +/- percent, -1 } }
			// or new double[][] { new double[] {param_num, (min_range, max_range) or (+/-
			// percent, -1) } ... {{-1 = default, (min_range, max_range) or (+/- percent,
			// -1)}}}
			// or new double[][] { new double[] {param_num, parameter_setting, ratio_adj_of}
			//
			// Fixed : null
			// Lookup : new double[][] {new double[] {top_entries,} ;

			param_range = new double[][] {       // 
				    new double[] { 3, 0.5, 1 }, // NG_trans U->R
					new double[] { 4, 0, 1, 3 }, // NG_trans R->U
					new double[] { 5, 0, 1, 3 }, // NG_trans U->O
					new double[] { 6, 0, 1, 5 }, // NG_trans O->U
					new double[] { 7, 0, 1, 3 }, // NG_trans R->O
					new double[] { 8, 0, 1, 7 }, // NG_trans O->R
					new double[] { 9, 0, 0.2 }, // NG_trans O->O
					new double[] { 10, 30, 200 }, // NG_dur U
					new double[] { 11, 30, 200 }, // NG_dur R
					new double[] { 12, 0, 1, 10 }, // NG_dur O
					new double[] { 15, 0.1, 0.2 }, // CT_trans U->R
					new double[] { 16, 0, 1, 15 }, // CT_trans R->U
					new double[] { 17, 0, 1, 15 }, // CT_trans U->O
					new double[] { 18, 230, 635 }, // CT_dur U
					new double[] { 19, 230, 635 }, // CT_dur R
					new double[] { 20, 0, 1, 18 }, // CT_dur O
					new double[] { 21, 0, 1 }, // Non-viability from treatment: NG_Rectal
					new double[] { 22, 0, 1 }, // Non-viability from treatment: NG_Orophraygneal
					new double[] { 23, 0, 1 }, // Non-viability from treatment: CT_Urethal
					new double[] { 24, 0, 1 }, // Non-viability from treatment: CT_Rectal
					new double[] { 25, 0, 1 }, // Non-viability from treatment: CT_Orophraygneal
					new double[] { -1, 0.10, -1 } // Default
			};
			numSimPerDir = 16;
			numDirPerMap = 5;
			break;
		case PROP_TYPE_SIM_BALI:
			outputDirName = "Bali_Baseline";
			idLine_str = "#PBS -N BalBas_%03d";
			templateDir = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Sim_Bali");
			sourcePBS = new File(
					"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Sim_Template_Bali.pbs");
			cMapDir = new File("C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\GenMap_BALI");
			lookup_file = null;

			dirFormat = "SimClusterModel_Transmission_Extra_%03d";
			seedHeader = "CMAP_SEED,SIM_SEED" + ",18_2048_32,18_16384_32,19_32_8,19_64_8,16_4_16,16_8_16" // NG
					+ ",18_131072_32,18_1048576_32,19_128_8,19_256_8,16_16_16,16_32_16" // CT
					+ ",16_1_16,19_31_8,18_4_32" // TP
					+ ",16_64_16" // HIV
					+ "";

			paramStr = ",170,360,0.80,0.02,0.968806,0.627861" // NG
					+ ",384,633,0.10,0.02,0.200000,0.120000" // CT
					+ ",0.11,0.050,450" // TP
					+ ",0.030" // HIV
					+ "";

			cMap_seedCopy_name = new String[0];
			direct_copy = new File[] { new File(templateDir, "simSpecificSim.prop"),
					new File(templateDir, "simSpecificSwitch.prop"), new File(templateDir, "dx_Bali.prop"), };
			// Usage: Pre-define range: new double[][] { new double[] { min_range, max_range
			// } } or new double[][] { new double[] { +/- percent, -1 } };
			// Fixed : null
			// Lookup : new double[][] {new double[] {top_entries,} ;
			param_range = new double[][] { new double[] { 0.10, -1 } };
			numMapInclMax = 125;
			numSimPerDir = 8;
			numDirPerMap = 1;

			break;
		case PROP_TYPE_SIM:
		default:
			extra_seed = 9;
			outputDirName = String.format("Baseline_Split_%d", extra_seed);
			idLine_str = String.format("#PBS -N SIM_S%d_%%03d", extra_seed);
			templateDir = new File(
					"C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Template\\SimClusterModel_Transmission_Sim");
			sourcePBS = new File(
					"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Sim_Template.pbs");
			cMapDir = new File("C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\Test_cMap");

			lookup_file = new File(
					"C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\SimClusterModel_Opt_Trend\\HPC Katana\\Current\\"
							+ "Opt_Shift\\CombinedOutput\\Opt_Summary.csv");

			dirFormat = "SimClusterModel_Transmission_Extra_%03d";
			seedHeader = "CMAP_SEED,SIM_SEED" + ",16_2_1_1,16_1_2_1,16_2_4_1" // P->V, V->P, P->R
					+ ",16_2_8_1,16_4_2_1,16_4_8_1" // P->O, R->P, R->O
					+ ",16_8_2_1,16_8_4_1,16_8_8_1" // O->P, O->R, O->O
					+ ",17_2_1,17_1_1,17_4_1,17_8_1" // Duration: penile, vaginal, rectal, orophargneal
					+ ",19_14_2" // Sym Seek MSM
					+ "";

			paramStr = "" + "0.818008,0.504899,0.346279" // P->V, V->P, P->R
					+ ",0.295040,0.093790,0.140796" // P->O, R->P, R->O
					+ ",0.007013,0.024413,0.004015" // O->P, O->R, O->O
					+ ",174.890750,235.528853,127.899082,78.581497" // Duration: penile, vaginal, rectal, phargneal
					+ ",0.523049" // Sym Seek MSM
					+ "";


			// Usage: Pre-define range: new double[][] { new double[] { min_range, max_range
			// } } or new double[][] { new double[] { +/- percent, -1 } };
			// Fixed : null
			// Lookup : new double[][] {new double[] {top_entries,} ;
			param_range = new double[][] { //
				    new double[] { 0, 0.6, 1 }, // P->V
					new double[] { 1, 0, 1, 0 }, // V->P
					new double[] { 2, 0, 0.9, 0 }, // P->R
					new double[] { 3, 0, 0.9, 2 }, // P->O
					new double[] { 4, 0, 1, 2 }, // R->P
					new double[] { 5, 0, 1, 2 }, // R->O
					new double[] { 6, 0, 1, 3 }, // O->P
					new double[] { 7, 0, 1, 3 }, // O->R
					new double[] { 8, 0, 0.01 }, // O->O
					new double[] { -1, 0.10, -1 } };
			numSimPerDir = 16;
			numDirPerMap = 7;
		}
		
		if(paramStr.startsWith(",")) {
			System.err.println("Warning: ParamStr start with ',' and possibly relies on old format - check validity of output.");
		}

		File targetDir = new File("C:\\Users\\bhui\\Documents\\Java_Test\\Prop_Gen");
		if (direct_copy == null) {
			direct_copy = new File[] { new File(templateDir, "simSpecificSim.prop"),
					new File(templateDir, "simSpecificSwitch.prop") };
		}

		MersenneTwisterRandomGenerator RNG = new MersenneTwisterRandomGenerator(rng_seed);

		for (int i = 0; i < extra_seed; i++) {
			RNG = new MersenneTwisterRandomGenerator(RNG.nextLong());
		}

		File[] cMaps = cMapDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isFile() && Pattern.matches(REGEX_CMAP, pathname.getName())
						&& !pathname.getName().startsWith("All_ContactMap_-1_");
			}
		});

		if (numMapInclMax < cMaps.length) {
			cMaps = Arrays.copyOf(cMaps, numMapInclMax);
		}

		ArrayList<File> generatedDirs = new ArrayList<>();

		int dir_counter = 0;

		boolean useParamLookup = lookup_file != null && param_range != null && param_range.length == 1
				&& param_range[0].length == 1;
		HashMap<String, double[]> param_lookup_map = null;
		RealDistribution[] param_dist = null;

		if (useParamLookup) {
			param_lookup_map = new HashMap<>();
			String[] def_str = paramStr.split(",");
			double[] def_val = new double[def_str.length - 2];
			long cMap_seed_ref = Long.parseLong(def_str[0]);
			long sim_seed_ref = Long.parseLong(def_str[1]);
			StringBuilder adj_paramStr = new StringBuilder();

			for (int p = 2; p < def_str.length; p++) {
				def_val[p - 2] = Double.parseDouble(def_str[p]);
				adj_paramStr.append(',');
				adj_paramStr.append(def_str[p]);
			}

			paramStr = adj_paramStr.toString();

			param_dist = new RealDistribution[def_val.length];

			String[] lines = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(lookup_file);
			double[] ref_line = null;

			int maxEntry = Math.min(lines.length, (int) param_range[0][0]);

			double[][] val_by_param = null;
			long[][] seeds = new long[2][maxEntry];

			for (int s = 0; s < maxEntry; s++) {
				String[] line_sp = lines[s].split(",");
				if (s == 0) {
					val_by_param = new double[def_val.length][maxEntry];
				} else {
					long cMap_seed = Long.parseLong(line_sp[1]);
					long sim_seed = Long.parseLong(line_sp[2]);
					seeds[0][s] = cMap_seed;
					seeds[1][s] = sim_seed;
					double[] param_val = new double[def_val.length];

					for (int p = 0; p < param_val.length; p++) {
						param_val[p] = Double.parseDouble(line_sp[p + 3]);
						val_by_param[p][s] = param_val[p];
					}
					if (cMap_seed == cMap_seed_ref && sim_seed == sim_seed_ref) {
						ref_line = param_val;
					}
				}
			}

			// Adjust parameter if ref_line is found, or use mean otherwise
			double[] meanVal = new double[def_val.length];
			for (int p = 0; p < meanVal.length; p++) {
				meanVal[p] = new DescriptiveStatistics(Arrays.copyOfRange(val_by_param[p], 1, val_by_param[p].length))
						.getMean();

			}
			for (int p = 0; p < val_by_param.length; p++) {
				for (int s = 0; s < val_by_param[p].length; s++) {
					if (ref_line != null) {
						val_by_param[p][s] *= def_val[p] / ref_line[p];
					} else {
						val_by_param[p][s] *= def_val[p] / meanVal[p];
					}
				}
			}

			// Look up and distribution

			for (int s = 1; s < maxEntry; s++) {
				String key = String.format("%d_%d", seeds[0][s], seeds[1][s]);
				if (!param_lookup_map.containsKey(key)) {
					double[] param_val = new double[def_val.length];
					for (int p = 0; p < param_val.length; p++) {
						param_val[p] = val_by_param[p][s];
					}
					param_lookup_map.put(key, param_val);
				}
			}

			for (int p = 0; p < param_dist.length; p++) {
				EmpiricalDistribution dist = new EmpiricalDistribution();
				dist.load(val_by_param[p]);
				param_dist[p] = dist;

			}

		}

		File genDirBase = new File(targetDir, outputDirName);
		genDirBase.mkdirs();

		PrintWriter pWri_Default_Setting = new PrintWriter(new File(genDirBase, "genDefaultSetting.txt"));
		pWri_Default_Setting.printf("Seed dirtecory generated at %1$td %1$tB, %1$tY %1$tH:%1$tM:%1$tS.\n", new Date());
		pWri_Default_Setting.printf("Output Dir = %s\n", outputDirName);
		pWri_Default_Setting.printf("DirFormat = %s\n", dirFormat);
		pWri_Default_Setting.printf("SeedHeader = %s\n", seedHeader);
		pWri_Default_Setting.printf("paramStr = %s\n", paramStr);
		pWri_Default_Setting.printf("RNG_Seed = %d\n", rng_seed);
		if (param_range != null) {
			pWri_Default_Setting.printf("paramStr = %s\n", Arrays.deepToString(param_range));
		}

		pWri_Default_Setting.close();

		for (File cMap : cMaps) {
			Matcher m = Pattern.compile(REGEX_CMAP).matcher(cMap.getName());
			m.matches();

			String cMapSeedStr = m.group(1);

			for (int i = 0; i < numDirPerMap; i++) {
				File genDir = new File(genDirBase, String.format(dirFormat, dir_counter));
				generatedDirs.add(genDir);
				if (genDir.exists()) {
					System.err.printf("Warning: Directory %s already existed in %s. Output directory NOT generated.\n",
							genDir.getName(), genDirBase.getAbsolutePath());
				} else {

					genDir.mkdirs();

					for (File copyFile : direct_copy) {
						if (copyFile.exists()) {
							File tarFile = new File(genDir, copyFile.getName());
							Files.copy(copyFile.toPath(), tarFile.toPath());
						}
					}

					for (String fName : cMap_seedCopy_name) {
						File srcFile = new File(cMapDir, String.format(fName, cMapSeedStr));
						if (!srcFile.isFile()) {
							System.err.printf("Warning! File %s not found.\n", srcFile.getAbsolutePath());
						} else {
							File tarFile = new File(genDir, srcFile.getName());
							Files.copy(srcFile.toPath(), tarFile.toPath());
						}

					}

					File seedFile = new File(genDir, "Seed_List.csv");
					PrintWriter pWri = new PrintWriter(seedFile);
					pWri.println(seedHeader);
					for (int j = 0; j < numSimPerDir; j++) {
						long sim_seed = RNG.nextLong();
						pWri.print(cMapSeedStr);
						pWri.print(',');
						pWri.print(sim_seed);
						if (param_range == null) {
							pWri.println(paramStr);
						} else if (useParamLookup) {
							String key = String.format("%s_%d", cMapSeedStr, sim_seed);
							if (param_lookup_map.containsKey(key)) {
								double[] ent = param_lookup_map.get(key);
								for (int p = 0; p < ent.length; p++) { // Param started at 3rd index
									pWri.print(',');
									pWri.print(ent[p]);
								}
							} else if (param_dist != null) {
								for (int p = 0; p < param_dist.length; p++) {
									pWri.print(',');
									pWri.print(param_dist[p].sample());

								}
							} else {
								// Use default
								pWri.print(paramStr);
							}
							pWri.println();

						} else {
							String[] param_str = paramStr.split(",");
							HashMap<Integer, Double> paramValMap = new HashMap<>();

							if (param_range.length == param_str.length || param_range.length == 1) {
								for (int r = 0; r < param_str.length; r++) {
									double paramVal;
									// Direct mapping or single mapping
									double[] range = r < param_range.length ? param_range[r] : param_range[0];
									if (range[1] < range[0]) {
										float param = Float.parseFloat(param_str[r + 1]);
										paramVal = param + param * 2 * range[0] * (RNG.nextFloat() - 0.5f);
									} else {
										paramVal = range[0] + RNG.nextFloat() * (range[1] - range[0]);
									}

									pWri.print(',');
									paramValMap.put(r, paramVal);
									pWri.print(paramVal);

								}
								pWri.println();
							} else {
								// Format: new double[][] { new double[] {param_num, parameter_setting,
								// ratio_adj_of}

								double[] paramVals = new double[param_str.length]; // NaN if need to be look up (i.e.
																					// based on other parameter)

								for (int r = 0; r < param_str.length; r++) {
									double[] range_r = null;
									for (double[] range : param_range) {
										if ((int) range[0] == r) {
											range_r = range;
											if (range_r.length > 3) {
												paramVals[r] = Double.NaN;
											}
										}
									}
									if (!Double.isNaN(paramVals[r])) {
										if (range_r == null) {
											range_r = param_range[param_range.length - 1];
										}
										if (range_r[1] > range_r[2]) {
											float defaultParam = Float.parseFloat(param_str[r]);
											paramVals[r] = defaultParam
													+ defaultParam * 2 * range_r[1] * (RNG.nextFloat() - 0.5f);
										} else {
											paramVals[r] = range_r[1] + RNG.nextFloat() * (range_r[2] - range_r[1]);
										}
										paramValMap.put(r, paramVals[r]);
									}
								}
								for (int r = 0; r < param_str.length; r++) {
									pWri.print(',');
									if (Double.isNaN(paramVals[r])) {
										// Look up
										for (double[] range : param_range) {
											if ((int) range[0] == r) {
												paramVals[r] = paramVals[(int) range[3]];
												paramVals[r] *= range[1] + RNG.nextDouble() * (range[2] - range[1]);
												break;
											}
										}
									}
									pWri.print(paramVals[r]);
								}
							}
							pWri.println();
						}
					}
					pWri.close();

					System.out.printf("Directory %s, utilising cMap %s, generated.\n", genDir.getName(),
							cMap.getName());
				}
				dir_counter++;

			}

		}
		// Generate PBS files
		File tarPBS = new File(targetDir, String.format("PBS_%s", outputDirName));
		tarPBS.mkdirs();

		String pattern_PBS_filename = "ClusterModel_Sim_%03d.pbs";
		String batch_script_filename = "batch_qsub";

		String[] line_pbs_arr = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(sourcePBS);

		System.out.printf("Number of line read from PBS = %d\n", line_pbs_arr.length);

		PrintWriter batchScriptWriter = new PrintWriter(new File(tarPBS, batch_script_filename)) {
			@Override
			public void println() {
				write('\n');
			}
		};

		batchScriptWriter.println("#!/bin/bash");
		batchScriptWriter.println("echo \"Batch submiting...\"");

		dir_counter = 0;

		for (File genDir : generatedDirs) {
			String pBS_FileName = String.format(pattern_PBS_filename, dir_counter);

			PrintWriter pbs_Writer = new PrintWriter(new File(tarPBS, pBS_FileName)) {
				@Override
				public void println() {
					write('\n');
				}
			};

			// Id line
			String id_line = String.format(idLine_str, dir_counter);
			String run_line = "";

			String optMethod, paramVal, paramLim, dirPath;

			switch (propType) {
			case PROP_TYPE_OPT_MSM:
				// Run line
				optMethod = "-opt_trend_fs";

				paramVal = "[" + "0.028000,180.0000,0.006750," // TP
						+ "0.549154,0.339017,0.231058,0.0666180,0.029145,0.012548,0.003501,61.141255,154.422144,108.275146,0.796314,0.443171," // NG
						+ "0.128574,0.049310,0.010380,384.95404,618.176299,289.239346" // CT
						+ "]";

				paramLim = "[[" + "0.005,30,0.001," // TP
						+ "0,0,0,0,0,0,0,30,30,30,0,0," // NG
						+ "0,0,0,14,14,14" // CT
						+ "],[" + "1,300,0.6," + "1,1,1,1,1,1,1,365,800,800,1,1," // NG
						+ "1,1,1,950,950,950" // CT
						+ "]]";

				dirPath = String.format("\"/srv/scratch/z2251912/%s/%s\"", outputDirName, genDir.getName());

				run_line = String.format("java -jar \"ClusterModel.jar\" %s %s %s %s" + " Seed_List.csv -nEval=1000",
						optMethod, dirPath, paramVal, paramLim);

				break;
			case PROP_TYPE_OPT_BALI:
				// Run line
				optMethod = "-opt_trend_fs";

				// All

				// Param 16 - Transmisson rate
				// Param 18 - Duration
				// Param 19 - Sym rate

				// NG :18_2048_32,18_16384_32,19_32_8,19_64_8,16_4_16,16_8_16,
				// CT :18_131072_32,18_1048576_32,19_128_8,19_256_8,16_16_16,16_32_16,
				// TP :18_4_32,19_31_8,16_1_16,
				// HIV:16_64_16

				paramVal = "[" + "167.761946,570.679869,0.950171,0.095177,0.93784,0.756719," // MG
						+ "47.117491,592.460657,0.352183,0.055024,0.967758,0.836822," // CT
						+ "360.000,0.002,0.0900," // TP
						+ "0.04451" // HIV
						+ "]";
				paramLim = "[[" + "14,60,0,0,0,0," // NG
						+ "14,60,0,0,0,0," // CT
						+ "21,0,0," // TP
						+ "0" // HIV
						+ "],[" + "180,600,1,1,1,1," // NG
						+ "600,600,1,1,1,1," // CT
						+ "600,1,1," // TP
						+ "1" // HIV
						+ "]]";

				dirPath = String.format("\"/srv/scratch/z2251912/%s/%s\"", outputDirName, genDir.getName());

				run_line = String.format("java -jar \"ClusterModel.jar\" %s %s %s %s" + " Seed_List.csv -nEval=1000",
						optMethod, dirPath, paramVal, paramLim);

				break;
			case PROP_TYPE_OPT_MSM_VIABILITY:
				// Run line
				optMethod = "-opt_trend_fs";

				paramVal = "["
						+ "0.549154,0.339017,0.231058,0.066618,0.029145,0.012548,0.003501,61.1412550,154.422144,108.275146," // NG
						+ "0.796314,0.443171,0.128574,0.049310,0.010380,384.954040,618.176299,289.239346," // CT
						+ "0.5,0.5,0.5,0.5" // Viability setting
						+ "]";
				paramLim = "[[0,0,0,0,0,0,0,30,30,30,0,0,0,0,0,30,30,30,0,0,0,0],"
						+ "[1,1,1,1,1,1,1,900,900,900,1,1,1,1,1,900,900,900,1,1,1,1]]";

				dirPath = String.format("\"/srv/scratch/z2251912/%s/%s\"", outputDirName, genDir.getName());

				run_line = String.format("java -jar \"ClusterModel.jar\" %s %s %s %s" + " Seed_List.csv -nEval=1000",
						optMethod, dirPath, paramVal, paramLim);

				break;
			case PROP_TYPE_OPT:
				// Run line
				optMethod = "-opt_trend_fs";

				paramVal = "[0.474899,0.818008,0.306279,0.292540,0.093790,0.140796,0.007013,0.024413,0.004015"
						+ ",174.890750,220.528853,147.899082,78.581497,0.523049]";
				paramLim = "[[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0," + "60,60,60,60,0"
						+ "],[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0," + "600,600,600,600,1]]";

				dirPath = String.format("\"/srv/scratch/z2251912/%s/%s\"", outputDirName, genDir.getName());

				run_line = String.format("java -jar \"ClusterModel.jar\" %s %s %s %s" + " Seed_List.csv -nEval=150",
						optMethod, dirPath, paramVal, paramLim);

				break;
			case PROP_TYPE_SIM_BALI:
			case PROP_TYPE_SIM_MSM_DOXY_PEP:
			case PROP_TYPE_SIM_MSM_VIABILITY:
			case PROP_TYPE_SIM:

				run_line = String.format(
						"java -jar \"ClusterModel.jar\" " + "-trans \"/srv/scratch/z2251912/All_Transmission/%s/%s\" "
								+ "-export_skip_backup -printProgress -seedMap=Seed_List.csv",
						outputDirName, genDir.getName());
				break;
			default:
				System.out.printf("PBS Generation for propType=%d be implemented. Exiting...\n", propType);
				System.exit(-1);
			}

			for (int n = 0; n < line_pbs_arr.length; n++) {
				String lineEnt = line_pbs_arr[n];
				switch (n) {
				case 5:
					lineEnt = id_line;
					break;
				case 9:
					lineEnt = String.format("#PBS -l ncpus=%d", numSimPerDir);
					break;
				case 14:
					lineEnt = String.format("#PBS -o ./%s.out", pBS_FileName);
					break;
				case 15:
					lineEnt = String.format("#PBS -e ./%s.err", pBS_FileName);
					break;
				case 23:
					lineEnt = run_line;
				default:

				}
				pbs_Writer.println(lineEnt);
			}

			pbs_Writer.close();

			batchScriptWriter.printf("qsub %s\n", pBS_FileName);
			dir_counter++;
		}

		batchScriptWriter.println("echo \"Submit completed!\"");
		batchScriptWriter.close();

		System.out.printf("PBS file generated from %d directories.\n", generatedDirs.size());

	}

}
