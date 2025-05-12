package sim;

import java.io.PrintStream;
import java.util.Properties;

import person.AbstractIndividualInterface;
import relationship.ContactMap;
import util.PropValUtils;

public abstract class Abstract_Runnable_ClusterModel_ContactMap_Generation extends Abstract_Runnable_ClusterModel {

	public static final Object[] DEFAULT_RUNNABLE_MAP_GEN_FIELDS = {
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE
			new int[] { 30, 30 + AbstractIndividualInterface.ONE_YEAR_INT },
			// RUNNABLE_FILED_EXPORT_FREQ - in milliseconds
			-1l, };

	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE = 0;
	public static final int RUNNABLE_FILED_EXPORT_FREQ = RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE + 1;
	public static final int LENGTH_RUNNABLE_MAP_GEN_FIELD = RUNNABLE_FILED_EXPORT_FREQ + 1;

	public static final String EXPORT_POP_FILENAME_PRRFIX = "EXPORT_POP_%d";
	public static final String EXPORT_POP_FILENAME = EXPORT_POP_FILENAME_PRRFIX + "_%d.obj";

	protected int numSnaps;
	protected int snap_dur;
	protected Object[] runnable_fields = new Object[LENGTH_RUNNABLE_MAP_GEN_FIELD];
	protected ContactMap[] gen_cMap = null;
	protected PrintStream[] printStatus = null;
	protected long mapSeed;

	protected boolean space_save = false;

	public Abstract_Runnable_ClusterModel_ContactMap_Generation(long mapSeed) {
		super();
		for (int i = 0; i < DEFAULT_RUNNABLE_MAP_GEN_FIELDS.length; i++) {
			runnable_fields[i] = DEFAULT_RUNNABLE_MAP_GEN_FIELDS[i];
		}
		this.mapSeed = mapSeed;

	}

	public boolean isSpace_save() {
		return space_save;
	}

	public void setSpace_save(boolean space_save) {
		this.space_save = space_save;
	}

	public long getMapSeed() {
		return mapSeed;
	}

	public void setNumSnaps(int numSnaps) {
		this.numSnaps = numSnaps;
	}

	public ContactMap[] getGen_cMap() {
		return gen_cMap;
	}

	public void setSnapFreq(int snapFreq) {
		this.snap_dur = snapFreq;
	}

	public void setPrintStatus(PrintStream[] printStatus) {
		this.printStatus = printStatus;
	}

	@Override
	public Object[] getRunnable_fields() {
		return runnable_fields;
	}

	public abstract void setRunnable_fields(Object[] simFields);

	public void loadRunnableFieldsFromProperties(Properties prop) {
		for (int i = 0; i < getRunnable_fields().length; i++) {
			String propName = String.format("%s%d", Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX, i);
			if (prop.containsKey(propName)) {
				String objStr = prop.getProperty(propName);
				Class<?> objClass = Object.class;

				if (getRunnable_fields()[i] != null) {
					objClass = getRunnable_fields()[i].getClass();
				} else {
					String str = String.format("%s%d", Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX_CLASS, i);
					if (prop.containsKey(str)) {
						str = prop.getProperty(str);
						try {
							objClass = Class.forName(str);
						} catch (ClassNotFoundException e) {
							e.printStackTrace();
						}
					}
				}
				getRunnable_fields()[i] = PropValUtils.propStrToObject(objStr, objClass);
			}

		}
	}

}
