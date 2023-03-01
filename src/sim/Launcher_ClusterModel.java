package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.InvalidPropertiesFormatException;

import optimisation.Optimisation_Factory;
import util.Util_Analyse_ClusterModel_Transmission_Output;
import util.Util_Analyse_ContactMap_Outputs;
import util.Util_Compare_ClusterModel_Transmission_Output;

public class Launcher_ClusterModel {

	public static void main(String[] args) throws InvalidPropertiesFormatException, IOException, InterruptedException {

		final String USAGE_INFO = String.format(
				"Usage: java %s <-gen, -trans, -opt, -analyse, -analyse_map or -compare> PROP_FILE_DIRECTORY <...>\n"
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
			} else if("-optGA".equals(flag)) {
				Optimisation_Factory.stable_prevalence_by_tranmission_fit_GA(Arrays.copyOfRange(args, 1, args.length));
			} else if ("-analyse".equals(flag)) {
				Util_Analyse_ClusterModel_Transmission_Output analysis = new Util_Analyse_ClusterModel_Transmission_Output();
				analysis.setBaseDir(new File(args[1]));
				analysis.analyse_outputs();
			} else if("-analyse_map".equals(flag)) {
				Util_Analyse_ContactMap_Outputs analysis = new Util_Analyse_ContactMap_Outputs(args[1], 
						Integer.parseInt(args[2]), Integer.parseInt(args[3]));
				analysis.analyse_output();
			} else if ("-compare".equals(flag)) {
				Util_Compare_ClusterModel_Transmission_Output.launch(Arrays.copyOfRange(args, 1, args.length));

			} else if ("-batch".equals(flag)) {
				File commands = new File(args[1]);
				BufferedReader reader = new BufferedReader(new FileReader(commands));
				String line;
				while ((line = reader.readLine()) != null) {
					String[] lines = line.split("\\s+(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", -1);
					for (int i = 0; i < lines.length; i++) {
						lines[i] = lines[i].replaceAll("\"", "");
					}
					Launcher_ClusterModel.main(lines);
				}
				reader.close();

			} else {
				System.out.println(USAGE_INFO);
				System.exit(0);
			}

		}

	}

}
