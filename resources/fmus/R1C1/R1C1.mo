within ;
model R1C1
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C=C, T(fixed=
          true, start=298.15))
    annotation (Placement(transformation(extent={{26,0},{46,20}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor thermalResistor(R=R)
              annotation (Placement(transformation(extent={{2,-10},{22,10}})));
  Modelica.Blocks.Interfaces.RealInput q "W"
    annotation (Placement(transformation(extent={{-126,20},{-86,60}})));
  Modelica.Blocks.Interfaces.RealInput Tout
    annotation (Placement(transformation(extent={{-126,-60},{-86,-20}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-38,-50},{-18,-30}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(extent={{-66,30},{-46,50}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{40,-50},{60,-30}})));
  Modelica.Blocks.Interfaces.RealOutput Tr
    annotation (Placement(transformation(extent={{98,-50},{118,-30}})));
  Modelica.Blocks.Interfaces.RealOutput qout "W"
    annotation (Placement(transformation(extent={{98,60},{118,80}})));
  parameter Modelica.SIunits.HeatCapacity C=1e7
    "Heat capacity of element (= cp*m)";
  parameter Modelica.SIunits.ThermalResistance R=0.01
    "Constant thermal resistance of material";
equation
  connect(heatCapacitor.port, thermalResistor.port_b)
    annotation (Line(points={{36,0},{22,0}}, color={191,0,0}));
  connect(prescribedTemperature.port, thermalResistor.port_a) annotation (Line(
        points={{-18,-40},{-8,-40},{-8,0},{2,0}}, color={191,0,0}));
  connect(q, prescribedHeatFlow.Q_flow)
    annotation (Line(points={{-106,40},{-66,40}}, color={0,0,127}));
  connect(prescribedHeatFlow.port, heatCapacitor.port) annotation (Line(points=
          {{-46,40},{50,40},{50,0},{36,0}}, color={191,0,0}));
  connect(heatCapacitor.port, temperatureSensor.port)
    annotation (Line(points={{36,0},{36,-40},{40,-40}}, color={191,0,0}));
  connect(Tout, prescribedTemperature.T)
    annotation (Line(points={{-106,-40},{-40,-40}}, color={0,0,127}));
  connect(temperatureSensor.T, Tr)
    annotation (Line(points={{60,-40},{108,-40}}, color={0,0,127}));
  connect(q, qout) annotation (Line(points={{-106,40},{-78,40},{-78,70},{108,70}},
        color={0,0,127}));
  annotation (uses(Modelica(version="3.2.2")));
end R1C1;
