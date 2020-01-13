package Web_Services;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.cloudbus.cloudsim.Host;
import org.cloudbus.cloudsim.Log;
import org.cloudbus.cloudsim.Vm;
import org.cloudbus.cloudsim.core.CloudSim;

public class Federate_VmAllocationPolicy extends org.cloudbus.cloudsim.VmAllocationPolicy {

	private final Map<String, Host> _vmTable = new HashMap<String, Host>();
	private final ListHost _hosts;

	public Federate_VmAllocationPolicy(List<? extends Host> list) {
		super(list);
		this._hosts = new ListHost(list);
	}

	@Override
	public boolean allocateHostForVm(Vm vm) {
		if (this._vmTable.containsKey(vm.getUid())) {
			return true;
		}

		boolean vm_allocated = false;

		Host host = this._hosts.next();
		if (host != null) {
			vm_allocated = this.allocateHostForVm(vm, host);
		}

		return vm_allocated;
	}

	@Override
	public boolean allocateHostForVm(Vm vm, Host host) {
		if (host != null && host.vmCreate(vm)) {
			_vmTable.put(vm.getUid(), host);
			Log.formatLine(
					"%.4f: VM #" + vm.getId() + " has been allocated to the host#" + host.getId() + " datacenter #"
							+ host.getDatacenter().getId() + "(" + host.getDatacenter().getName() + ") #",
					CloudSim.clock());
			return true;
		}
		return false;
	}

	@Override
	public List<Map<String, Object>> optimizeAllocation(List<? extends Vm> vmList) {
		return null;
	}

	@Override
	public void deallocateHostForVm(Vm vm) {
		Host host = this._vmTable.remove(vm.getUid());

		if (host != null) {
			host.vmDestroy(vm);
		}
	}

	@Override
	public Host getHost(Vm vm) {
		return this._vmTable.get(vm.getUid());
	}

	@Override
	public Host getHost(int vmId, int userId) {
		return this._vmTable.get(Vm.getUid(userId, vmId));
	}
}
