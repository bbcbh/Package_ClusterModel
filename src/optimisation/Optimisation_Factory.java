package optimisation;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import population.Population_Bridging;
import sim.Abstract_Runnable_ClusterModel;
import sim.Abstract_Runnable_ClusterModel_ContactMap_Generation;
import sim.Simulation_ClusterModelGeneration;
import sim.Simulation_ClusterModelTransmission;

public class Optimisation_Factory {

	public static final String OPT_OUTPUT_PREFIX_CMAP = "CMAP    = ";
	public static final String OPT_OUTPUT_PREFIX_SIMSEED = "SimSeed = ";
	public static final String OPT_OUTPUT_PREFIX_PARAM = "Param   = ";
	public static final String OPT_OUTPUT_PREFIX_RESIDUE = "Residue = ";
	public static final int RUNNABLE_OFFSET = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
				+ Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
				+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD;

	public static void recursiveRunnableFieldReplace(Object runnableField, int param_index, double[] param_val_all,
			String[] param_setting_all, int setting_level) {
		BigInteger arraySel = new BigInteger(param_setting_all[setting_level]);
		double offset = 0;

		if (runnableField instanceof int[] || runnableField instanceof float[] || runnableField instanceof double[]) {
			Matcher m = OptFittingFunction.POP_PROP_OPT_PARAM_FIT_SETTING_DIFF_FORMAT
					.matcher(param_setting_all[param_setting_all.length - 1]);
			if (m.find()) {
				int offsetIndex = Integer.parseInt((m.group(1)));
				offset = param_val_all[offsetIndex];
			}

			if (runnableField instanceof int[]) {
				int[] val_int_array = (int[]) runnableField;
				for (int i = 0; i < val_int_array.length; i++) {
					if (arraySel.testBit(i)) {
						val_int_array[i] = (int) Math.round(offset + param_val_all[param_index]);
					}
				}
			} else if (runnableField instanceof float[]) {
				float[] val_float_array = (float[]) runnableField;
				for (int i = 0; i < val_float_array.length; i++) {
					if (arraySel.testBit(i)) {
						val_float_array[i] = (float) (offset + param_val_all[param_index]);
					}
				}
			} else if (runnableField instanceof double[]) {
				double[] val_double_array = (double[]) runnableField;
				for (int i = 0; i < val_double_array.length; i++) {
					if (arraySel.testBit(i)) {
						val_double_array[i] = offset + param_val_all[param_index];
					}
				}
			}

		} else {
			if (runnableField.getClass().isArray()) {
				Object[] obj_array = (Object[]) runnableField;
				for (int i = 0; i < obj_array.length; i++) {
					if (arraySel.testBit(i)) {
						recursiveRunnableFieldReplace(obj_array[i], param_index, param_val_all, param_setting_all,
								setting_level + 1);

					}
				}
			} else {
				System.err.printf("Class contructor for %s not supported (wrong param number?). Exiting.\n",
						runnableField.getClass().getName());
				System.exit(1);
			}

		}
	}

	public static void setOptParamInRunnable_Transfrom(Abstract_Runnable_ClusterModel target_runnable, String transform_str, HashMap<String, Double> param_map, ArrayList<Integer> field_to_update) {
		// Fill transfrom_ents
		ArrayList<String[]> transfrom_ents = new ArrayList<>();
		Pattern pattern_braceEntry = Pattern.compile("\\[([^\\[\\]]*)\\]");
	
		Matcher m = pattern_braceEntry.matcher(transform_str.substring(1, transform_str.length() - 1));
		while (m.find()) {
			String ent = m.group(1);
			transfrom_ents.add(ent.split(","));
		}
	
		for (String[] transfrom_ent : transfrom_ents) {
			String transfrom_param_name = transfrom_ent[0];
			String[] param_setting_arr = transfrom_param_name.split("_");
			int param_name_index = Integer.parseInt(param_setting_arr[0]);
			int field_id = param_name_index - RUNNABLE_OFFSET;
			Object val = target_runnable.getRunnable_fields()[field_id];
	
			if (val != null) {
				double transformed_val = Double.parseDouble(transfrom_ent[transfrom_ent.length - 1]);
				int pt = 1;
				while (pt < transfrom_ent.length - 1) {
					String srcParam = transfrom_ent[pt];
					double base = 1;
					if (srcParam.startsWith("*")) {
						if (param_map.containsKey(transfrom_param_name)) {
							base = param_map.get(transfrom_param_name).doubleValue();
						} else {
							Object baseVal = val;
							for (int i = 1; i < param_setting_arr.length; i++) {
								int incIndex = Integer.parseInt(param_setting_arr[i]);
								if (incIndex != 0) {
									int shiftPt = 0;
									while ((incIndex & 1 << shiftPt) == 0) {
										shiftPt++;
									}
									if (baseVal instanceof Object[]) {
										baseVal = ((Object[]) baseVal)[shiftPt];
									} else {
										base = ((double[]) baseVal)[shiftPt];
									}
								}
							}
							param_map.put(transfrom_param_name, base);
						}
						srcParam = srcParam.substring(1);
					}
	
					double paramVal;
					if (param_map.containsKey(srcParam)) {
						paramVal = param_map.get(srcParam).doubleValue();
					} else {
						paramVal = Double.NaN;
						String[] srcParamSplit = srcParam.split("_");
						Object srcVal = target_runnable.getRunnable_fields()[Integer.parseInt(srcParamSplit[0])
								- RUNNABLE_OFFSET];
	
						for (int i = 1; i < srcParamSplit.length; i++) {
							int incIndex = Integer.parseInt(srcParamSplit[i]);
							if (incIndex != 0) {
								int shiftPt = 0;
								while ((incIndex & 1 << shiftPt) == 0) {
									shiftPt++;
								}
								if (srcVal instanceof Object[]) {
									srcVal = ((Object[]) srcVal)[shiftPt];
								} else if (srcVal instanceof float[]) {
									paramVal = ((float[]) srcVal)[shiftPt];
								} else {
									paramVal = ((double[]) srcVal)[shiftPt];
								}
							}
						}
						param_map.put(srcParam, paramVal);
					}
	
					transformed_val += base * paramVal * Double.parseDouble(transfrom_ent[pt + 1]);
	
					pt += 2;
				}
				int setting_level = 1;
				recursiveRunnableFieldReplace(val, 0, new double[] { transformed_val }, param_setting_arr,
						setting_level);
				pt = Collections.binarySearch(field_to_update, field_id);
				if (pt < 0) {
					field_to_update.add(~pt, field_id);
				}
				param_map.put(transfrom_param_name, transformed_val);
			}
		}
	}

	public Optimisation_Factory() {
		super();
	}

}