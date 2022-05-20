package sim;

import java.io.IOException;
import java.util.Arrays;
import java.util.InvalidPropertiesFormatException;

public class Launcher_ClusterModel {

	public static void main(String[] args) throws InvalidPropertiesFormatException, IOException, InterruptedException {
		
		final String USAGE_INFO = String.format("Usage: java %s <-gen or -trans> PROP_FILE_DIRECTORY <...>",
				Launcher_ClusterModel.class.getName());
		
		if (args.length < 1) {
			System.out.println(USAGE_INFO);
			System.exit(0);
		} else {
			String flag = args[0];

			if ("-gen".equals(flag)) {
				Simulation_ClusterModelGeneration.launch(Arrays.copyOfRange(args, 1, args.length));
			} else if ("-trans".equals(flag)) {
				Simulation_ClusterModelTransmission.launch(Arrays.copyOfRange(args, 1, args.length));
			} else {
				System.out.println(USAGE_INFO);
				System.exit(0);
			}

		}

	}

}
