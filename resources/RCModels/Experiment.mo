within RCModels;
model Experiment
  R1C1 r1C1_1 annotation (Placement(transformation(extent={{-12,-8},{8,12}})));
  Modelica.Blocks.Sources.Constant const(k=1000)
    annotation (Placement(transformation(extent={{-72,-6},{-52,14}})));
  Modelica.Blocks.Sources.Constant const1(k=1)
    annotation (Placement(transformation(extent={{-72,-40},{-52,-20}})));
  Modelica.Blocks.Sources.Constant const2(k=273.15)
    annotation (Placement(transformation(extent={{-72,26},{-52,46}})));
equation
  connect(const.y, r1C1_1.qin) annotation (Line(points={{-51,4},{-32,4},{-32,-2},
          {-12.2,-2}}, color={0,0,127}));
  connect(const1.y, r1C1_1.price) annotation (Line(points={{-51,-30},{-32,-30},
          {-32,-5.4},{-12.2,-5.4}}, color={0,0,127}));
  connect(const2.y, r1C1_1.Tamb) annotation (Line(points={{-51,36},{-26,36},{
          -26,3.2},{-12.2,3.2}}, color={0,0,127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Experiment;
