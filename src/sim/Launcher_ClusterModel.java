package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.InvalidPropertiesFormatException;
import java.util.regex.Pattern;

import optimisation.Optimisation_Factory;
import util.Util_Analyse_ClusterModel_Transmission_Output;
import util.Util_Analyse_ContactMap_Outputs;
import util.Util_Combine_ContactMap;
import util.Util_Compare_ClusterModel_Transmission_Output;

public class Launcher_ClusterModel {

	public static void main(String[] args) throws InvalidPropertiesFormatException, IOException, InterruptedException {

		final String USAGE_INFO = String.format(
				"Usage: java %s <-gen, -trans, -opt, -analyse, -analyse_rx -analyse_map, -combine_map, -compare -clean_up> PROP_FILE_DIRECTORY <...>\n"
						+ "or    java %s <-batch> COMMAND_AS_TEXT",
				Launcher_ClusterModel.class.getName(), Launcher_ClusterModel.class.getName());

		if (args.length < 1) {
			System.out.println(USAGE_INFO);
			System.exit(0);
		} else {
			String flag = args[0];

			if ("-gen".equals(flag)) {
				Simulation_ClusterModelGeneration.launch(Arrays.copyOfRange(args, 1, args.length));
			} else if ("-trans".equals(flag)) {
				Simulation_ClusterModelTransmission.launch(Arrays.copyOfRange(args, 1, args.length));
			} else if ("-opt".equals(flag)) {
				Optimisation_Factory
						.stable_prevalence_by_tranmission_fit_Simplex(Arrays.copyOfRange(args, 1, args.length));
			} else if ("-optGA".equals(flag)) {
				Optimisation_Factory.stable_prevalence_by_tranmission_fit_GA(Arrays.copyOfRange(args, 1, args.length));
			} else if ("-analyse".equals(flag)) {
				Util_Analyse_ClusterModel_Transmission_Output analysis = new Util_Analyse_ClusterModel_Transmission_Output();
				File dir = new File(args[1]);
				analysis.setBaseDir(dir);				
				System.out.printf("=== %s ===\n", dir.getName());
				if (args.length > 2) {
					analysis.setSkipAnalysis(Integer.parseInt(args[2]));
				}
				analysis.analyse_outputs();
			} else if ("-analyse_rx".equals(flag)) {
				File baseDir = new File(args[1]);
				if (args.length < 2) {
					System.out.printf("Usage: java %s -analyse_rx FILE_DIRECTORY PATTERN\n" );
					System.exit(0);
				} else {
					Pattern dirPattern = Pattern.compile(args[2]);
					File[] candidateDir = baseDir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return pathname.isDirectory() && dirPattern.matcher(pathname.getName()).matches();
						}
					});
					System.out.printf("Number of diretories to be analyse = %d\n", candidateDir.length);
					for (File dir : candidateDir) {
						System.out.printf("=== %s ===\n", dir.getName());
						Util_Analyse_ClusterModel_Transmission_Output analysis = new Util_Analyse_ClusterModel_Transmission_Output();
						analysis.setBaseDir(dir);
						if (args.length > 3) {
							analysis.setSkipAnalysis(Integer.parseInt(args[3]));
						}
						analysis.analyse_outputs();												
						Util_Analyse_ClusterModel_Transmission_Output.cleanUpOutputDir(dir);

					}
					System.out.printf("%d directories analysed and cleaned up.\n", candidateDir.length);
					
				}

			} else if ("-analyse_map".equals(flag)) {
				Util_Analyse_ContactMap_Outputs analysis = new Util_Analyse_ContactMap_Outputs(args[1],
						Integer.parseInt(args[2]), Integer.parseInt(args[3]));
				analysis.analyse_output();
			} else if ("-combine_map".equals(flag)) {
				Util_Combine_ContactMap combine = new Util_Combine_ContactMap(args[1], args[2], args[3]);
				combine.combineMaps();
			} else if ("-compare".equals(flag)) {
				Util_Compare_ClusterModel_Transmission_Output.launch(Arrays.copyOfRange(args, 1, args.length));
			}else if ("-clean_up".equals(flag)) {
				Util_Analyse_ClusterModel_Transmission_Output.cleanUpOutputDir(new File(args[1]));
			} else if ("-batch".equals(flag)) {
				File commands = new File(args[1]);
				BufferedReader reader = new BufferedReader(new FileReader(commands));
				String line;
				while ((line = reader.readLine()) != null) {
					if (!line.startsWith("//")) {
						String[] lines = line.split("\\s+(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", -1);
						for (int i = 0; i < lines.length; i++) {
							lines[i] = lines[i].replaceAll("\"", "");
						}
						Launcher_ClusterModel.main(lines);
					}
				}
				reader.close();

			} else {
				System.out.println(USAGE_INFO);
				System.exit(0);
			}

		}

	}

}
