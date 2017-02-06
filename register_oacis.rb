repo_dir = File.expand_path(File.dirname(__FILE__))

localhost = Host.find_by_name("localhost")

sim_params = {
  name: "DynamicalGraphModel",
  command: "#{repo_dir}/run.sh",
  support_input_json: false,
  print_version_command: "cd #{repo_dir} && git describe --always",
  parameter_definitions: [
    {key: "c", type: "Float", default: 0.2, description: "connection probability"},
    {key: "t_init", type: "Integer", default: 1024, description: "initial warm up duration"},
    {key: "t_measure", type: "Integer", default: 65536, description: "simulation duration"}
  ],
  description: "dynamical graph model",
  executable_on: [ localhost ]
}

sim_name = sim_params[:name]
if Simulator.where(name: sim_name).exists?
  puts "simulator #{sim_name} already exists" 
else
  sim = Simulator.create!(sim_params)
end

