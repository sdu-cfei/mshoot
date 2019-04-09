within ;
model R2C2
  Modelica.Blocks.Interfaces.RealInput u2
    annotation (Placement(transformation(extent={{-126,20},{-86,60}})));
  Modelica.Blocks.Interfaces.RealInput u1
    annotation (Placement(transformation(extent={{-126,-40},{-86,0}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor2(C=C2)
    annotation (Placement(transformation(extent={{20,40},{40,60}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor3(C=C3)
    annotation (Placement(transformation(extent={{70,40},{90,60}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor thermalResistor2(R=
        R2) annotation (Placement(transformation(extent={{-4,30},{16,50}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor thermalResistor3(R=
        R3) annotation (Placement(transformation(extent={{46,30},{66,50}})));
  parameter Modelica.SIunits.HeatCapacity C1=750000
    "Heat capacity of element (= cp*m)";
  parameter Modelica.SIunits.HeatCapacity C2=100000
    "Heat capacity of element (= cp*m)";
  parameter Modelica.SIunits.HeatCapacity C3=50000
    "Heat capacity of element (= cp*m)";

  parameter Modelica.SIunits.ThermalResistance R1=0.01
    "Constant thermal resistance of material";
  parameter Modelica.SIunits.ThermalResistance R2=0.01
    "Constant thermal resistance of material";
  parameter Modelica.SIunits.ThermalResistance R3=0.01
    "Constant thermal resistance of material";

  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{-78,30},{-58,50}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
    annotation (Placement(transformation(extent={{-78,-30},{-58,-10}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={80,16})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{96,-30},{116,-10}})));

  Modelica.Blocks.Interfaces.RealOutput u1_y
    annotation (Placement(transformation(extent={{96,-58},{116,-38}})));
equation
  connect(heatCapacitor2.port,thermalResistor2. port_b)
    annotation (Line(points={{30,40},{16,40}}, color={191,0,0}));
  connect(heatCapacitor2.port,thermalResistor3. port_a)
    annotation (Line(points={{30,40},{46,40}}, color={191,0,0}));
  connect(heatCapacitor3.port,thermalResistor3. port_b)
    annotation (Line(points={{80,40},{66,40}}, color={191,0,0}));
  connect(u2, prescribedHeatFlow.Q_flow)
    annotation (Line(points={{-106,40},{-78,40}}, color={0,0,127}));
  connect(u1, prescribedHeatFlow1.Q_flow)
    annotation (Line(points={{-106,-20},{-78,-20}}, color={0,0,127}));
  connect(prescribedHeatFlow1.port,heatCapacitor2. port)
    annotation (Line(points={{-58,-20},{30,-20},{30,40}}, color={191,0,0}));
  connect(heatCapacitor3.port, temperatureSensor.port)
    annotation (Line(points={{80,40},{80,26}}, color={191,0,0}));
  connect(temperatureSensor.T, y)
    annotation (Line(points={{80,6},{80,-20},{106,-20}}, color={0,0,127}));
  connect(y, y)
    annotation (Line(points={{106,-20},{106,-20}}, color={0,0,127}));
  connect(prescribedHeatFlow.port, thermalResistor2.port_a)
    annotation (Line(points={{-58,40},{-4,40}}, color={191,0,0}));
  connect(u1, u1_y) annotation (Line(points={{-106,-20},{-82,-20},{-82,-48},{
          106,-48}}, color={0,0,127}));
  annotation (uses(Modelica(version="3.2.2")));
end R2C2;
