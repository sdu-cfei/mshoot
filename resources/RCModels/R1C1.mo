within RCModels;
model R1C1
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitorC(C=C, T(
        fixed=true, start=293.15))
    annotation (Placement(transformation(extent={{20,12},{40,32}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor thermalResistorExt(R=R)
    annotation (Placement(transformation(extent={{-16,2},{4,22}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-60,2},{-40,22}})));
  Modelica.Blocks.Interfaces.RealInput Tamb
    annotation (Placement(transformation(extent={{-122,-8},{-82,32}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{-60,-50},{-40,-30}})));
  Modelica.Blocks.Interfaces.RealInput qin
    annotation (Placement(transformation(extent={{-122,-60},{-82,-20}})));
  parameter Modelica.SIunits.HeatCapacity C=1e6
    "Heat capacity of element (= cp*m)";
  parameter Modelica.SIunits.ThermalResistance R=0.01
    "Constant thermal resistance of material";
  Modelica.Blocks.Interfaces.RealOutput T
    annotation (Placement(transformation(extent={{96,2},{116,22}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{52,2},{72,22}})));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
    annotation (Placement(transformation(extent={{-6,-50},{14,-30}})));
  Modelica.Blocks.Interfaces.RealOutput q
    annotation (Placement(transformation(extent={{96,-60},{116,-40}})));
  Modelica.Blocks.Interfaces.RealInput price
    "\"Price per kWh of energy per hour\""
    annotation (Placement(transformation(extent={{-122,-94},{-82,-54}})));
  Modelica.Blocks.Math.Product product
    annotation (Placement(transformation(extent={{0,-84},{20,-64}})));
  Modelica.Blocks.Continuous.Integrator integrator(k=1/3600/1000)
    annotation (Placement(transformation(extent={{42,-84},{62,-64}})));
  Modelica.Blocks.Interfaces.RealOutput p
    annotation (Placement(transformation(extent={{96,-84},{116,-64}})));
equation
  connect(heatCapacitorC.port, thermalResistorExt.port_b)
    annotation (Line(points={{30,12},{4,12}},
                                            color={191,0,0}));
  connect(thermalResistorExt.port_a, prescribedTemperature.port)
    annotation (Line(points={{-16,12},{-40,12}},
                                               color={191,0,0}));
  connect(prescribedTemperature.T, Tamb)
    annotation (Line(points={{-62,12},{-102,12}},
                                                color={0,0,127}));
  connect(prescribedHeatFlow.Q_flow, qin)
    annotation (Line(points={{-60,-40},{-102,-40}}, color={0,0,127}));
  connect(heatCapacitorC.port, temperatureSensor.port)
    annotation (Line(points={{30,12},{52,12}},
                                             color={191,0,0}));
  connect(T, temperatureSensor.T)
    annotation (Line(points={{106,12},{72,12}},
                                              color={0,0,127}));
  connect(prescribedHeatFlow.port, heatFlowSensor.port_a)
    annotation (Line(points={{-40,-40},{-6,-40}}, color={191,0,0}));
  connect(heatFlowSensor.port_b, heatCapacitorC.port)
    annotation (Line(points={{14,-40},{30,-40},{30,12}},color={191,0,0}));
  connect(heatFlowSensor.Q_flow, q)
    annotation (Line(points={{4,-50},{106,-50}},         color={0,0,127}));
  connect(heatFlowSensor.Q_flow, product.u1) annotation (Line(points={{4,-50},{
          -14,-50},{-14,-68},{-2,-68}}, color={0,0,127}));
  connect(price, product.u2) annotation (Line(points={{-102,-74},{-54,-74},{-54,
          -80},{-2,-80}}, color={0,0,127}));
  connect(product.y, integrator.u)
    annotation (Line(points={{21,-74},{40,-74}}, color={0,0,127}));
  connect(integrator.y, p)
    annotation (Line(points={{63,-74},{106,-74}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)));
end R1C1;
