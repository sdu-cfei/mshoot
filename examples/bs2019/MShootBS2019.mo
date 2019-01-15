within ;
package MShootBS2019
  model ZoneR1C1 "Simple zone"

    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real maxHeat=5000 "Heating power of radiators [W]";
    parameter Real Tve=21 "Ventilation air temperature";
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput T "[K]"
      annotation (Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-120,80},{-100,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{62,110},{82,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Blocks.Math.Product hmltp
      annotation (Placement(transformation(extent={{-62,18},{-42,38}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-14,18},{6,38}})));
    Modelica.Blocks.Tables.CombiTable1D MetabolicHeat(table=[293.15,84.; 325.15,0.])
      annotation (Placement(transformation(extent={{-156,-2},{-132,22}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-146},{-4,-126}})));
    Modelica.Blocks.Interfaces.RealOutput Qr "Heat supplied by radiators [Wh]"
      annotation (Placement(transformation(extent={{216,-170},{236,-150}}),
          iconTransformation(extent={{214,-132},{234,-112}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
      annotation (Placement(transformation(extent={{-166,-54},{-146,-34}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-204,80},{-184,100}})));
    Modelica.Blocks.Math.Gain occeffectiv(k=occheff)
      annotation (Placement(transformation(extent={{-116,0},{-96,20}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cair(
                                    C=tmass*1.2*1005*Vi, T(fixed=false))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    AirMix                         airMix
      annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-170},{170,-150}})));
    Modelica.Blocks.Interfaces.RealInput occ "Number of occupants"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Interfaces.RealOutput vetot "Total airflow supply [m3]"
      annotation (Placement(transformation(extent={{216,-88},{236,-68}}),
          iconTransformation(extent={{214,-64},{234,-44}})));
    Modelica.Blocks.Continuous.Integrator integrator1(k=1)
      annotation (Placement(transformation(extent={{150,-88},{170,-68}})));
    Modelica.Blocks.Math.Max max1
      annotation (Placement(transformation(extent={{-130,-48},{-110,-28}})));
    Modelica.Blocks.Math.Gain ventilation(k=maxVent)
      annotation (Placement(transformation(extent={{-88,-88},{-68,-68}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-88,-146},{-68,-126}})));
    Modelica.Blocks.Interfaces.RealOutput ve "Ventilation airflow rate [m3/h]"
      annotation (Placement(transformation(extent={{216,-128},{236,-108}}),
          iconTransformation(extent={{214,-90},{234,-70}})));
    Modelica.Blocks.Interfaces.RealOutput qr "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-210},{236,-190}}),
          iconTransformation(extent={{214,-160},{234,-140}})));
    Modelica.Blocks.Math.Gain scale(k=0.01) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-150,-78})));
    Modelica.Blocks.Interfaces.RealInput dpos "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-242,-92},{-214,-64}}),
          iconTransformation(extent={{-234,-60},{-206,-32}})));
    Modelica.Blocks.Sources.Constant const2(k=Tve)
      annotation (Placement(transformation(extent={{-202,-54},{-182,-34}})));

    Modelica.Blocks.Interfaces.RealInput vpos "Radiator valve position [%]"
      annotation (Placement(transformation(extent={{-240,-150},{-212,-122}}),
          iconTransformation(extent={{-232,-152},{-204,-124}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{16,-60},{36,-40}})));
    Modelica.Blocks.Interfaces.RealOutput qv "Ventilation heat load [W]"
      annotation (Placement(transformation(extent={{216,24},{236,44}}),
          iconTransformation(extent={{214,2},{234,22}})));
    Modelica.Blocks.Continuous.Integrator integrator2(k=1)
      annotation (Placement(transformation(extent={{152,-24},{172,-4}})));
    Modelica.Blocks.Interfaces.RealOutput Qv "Ventilation heat load [Wh]"
      annotation (Placement(transformation(extent={{216,-24},{236,-4}}),
          iconTransformation(extent={{214,-26},{234,-6}})));
    Modelica.Blocks.Math.Gain scale1(k=0.01)
                                            annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-150,-136})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(hmltp.y, occHeatGain.Q_flow)
      annotation (Line(points={{-41,28},{-14,28}},          color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-24,-136}},                        color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-120,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{82,120},{180,120}},            color={0,0,127}));
    connect(T, fromKelvin.Celsius)
      annotation (Line(points={{224,120},{210,120},{203,120}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-183,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-206,90},{-226,90}},           color={0,0,127}));
    connect(hmltp.u2, occeffectiv.y) annotation (Line(points={{-64,22},{-80,22},{-80,
            10},{-95,10}},          color={0,0,127}));
    connect(MetabolicHeat.y[1], occeffectiv.u) annotation (Line(points={{-130.8,10},
            {-118,10}},                   color={0,0,127}));
    connect(temperatureSensor.T, MetabolicHeat.u[1]) annotation (Line(points={{82,120},
            {90,120},{90,-14},{-168,-14},{-168,10},{-158.4,10}},
          color={0,0,127}));
    connect(re.port_b, cair.port) annotation (Line(points={{-100,90},{-80,90},{
            -80,120},{46,120}}, color={191,0,0}));
    connect(prescribedHeatFlow.port, cair.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.port, cair.port)
      annotation (Line(points={{62,120},{62,120},{46,120}}, color={191,0,0}));
    connect(occHeatGain.port, cair.port)
      annotation (Line(points={{6,28},{46,28},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{82,120},{
            90,120},{90,-66},{0,-66},{0,-60.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, cair.port) annotation (Line(points={{-4,-136},
            {46,-136},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,-146},
            {-14,-160},{148,-160}},                      color={0,0,127}));
    connect(Qr, integrator.y)
      annotation (Line(points={{226,-160},{171,-160}}, color={0,0,127}));
    connect(occ, hmltp.u1) annotation (Line(points={{-226,34},{-64,34}},
                  color={0,0,127}));
    connect(integrator1.y, vetot)
      annotation (Line(points={{171,-78},{176,-78},{226,-78}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, max1.u1) annotation (Line(points={{-183,90},{-174,
            90},{-174,-32},{-132,-32}}, color={0,0,127}));
    connect(max1.y, airMix.Tve) annotation (Line(points={{-109,-38},{-40,-38},{
            -40,-46.6},{-10.8,-46.6}}, color={0,0,127}));
    connect(ventilation.y, airMix.Vve) annotation (Line(points={{-67,-78},{-28,-78},
            {-28,-55},{-11,-55}},      color={0,0,127}));
    connect(ventilation.y, integrator1.u)
      annotation (Line(points={{-67,-78},{-67,-78},{148,-78}},
                                                     color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,-136},
            {-67,-136}},                         color={0,0,127}));
    connect(ventilation.y, ve) annotation (Line(points={{-67,-78},{126.5,-78},{
            126.5,-118},{226,-118}}, color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qr) annotation (Line(points={{-14,-146},{-14,
            -200},{226,-200}}, color={0,0,127}));
    connect(ventilation.u, scale.y)
      annotation (Line(points={{-90,-78},{-139,-78}}, color={0,0,127}));
    connect(dpos, scale.u)
      annotation (Line(points={{-228,-78},{-162,-78}}, color={0,0,127}));
    connect(max1.u2, toKelvin.Kelvin)
      annotation (Line(points={{-132,-44},{-145,-44}}, color={0,0,127}));
    connect(toKelvin.Celsius, const2.y)
      annotation (Line(points={{-168,-44},{-181,-44}}, color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{10,-50},{16,-50}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, cair.port)
      annotation (Line(points={{36,-50},{46,-50},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qv) annotation (Line(points={{26,-60},{26,
            -72},{120,-72},{120,34},{226,34}}, color={0,0,127}));
    connect(integrator2.y, Qv)
      annotation (Line(points={{173,-14},{226,-14}}, color={0,0,127}));
    connect(integrator2.u, qv) annotation (Line(points={{150,-14},{120,-14},{
            120,34},{226,34}}, color={0,0,127}));
    connect(vpos, scale1.u)
      annotation (Line(points={{-226,-136},{-162,-136}}, color={0,0,127}));
    connect(scale1.y, heating.u)
      annotation (Line(points={{-139,-136},{-90,-136}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1)),
      experiment(
        Tolerance=1e-09,
        __Dymola_fixedstepsize=1,
        __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end ZoneR1C1;

  model ZoneR1C1_compact "Simple zone"

    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real maxHeat=5000 "Heating power of radiators [W]";
    parameter Real Tve=21 "Ventilation air temperature";
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-180,140},{-152,168}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-180,96},{-152,124}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput T "[K]"
      annotation (Placement(transformation(extent={{160,130},{180,150}}),
          iconTransformation(extent={{160,130},{180,150}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-60,100},{-40,120}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{62,130},{82,150}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-136,144},{-116,164}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-98,144},{-78,164}})));
    Modelica.Blocks.Math.Product hmltp
      annotation (Placement(transformation(extent={{-26,54},{-6,74}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{6,54},{26,74}})));
    Modelica.Blocks.Tables.CombiTable1D MetabolicHeat(table=[293.15,84.; 325.15,0.])
      annotation (Placement(transformation(extent={{-98,34},{-74,58}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-100,100},{-80,120}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{4,-110},{24,-90}})));
    Modelica.Blocks.Interfaces.RealOutput Qr "Heat supplied by radiators [Wh]"
      annotation (Placement(transformation(extent={{160,-130},{180,-110}}),
          iconTransformation(extent={{214,-132},{234,-112}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-34,-110},{-14,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{104,130},{124,150}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
      annotation (Placement(transformation(extent={{-118,-18},{-98,2}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-144,100},{-124,120}})));
    Modelica.Blocks.Math.Gain occeffectiv(k=occheff)
      annotation (Placement(transformation(extent={{-64,36},{-44,56}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cair(
                                    C=tmass*1.2*1005*Vi, T(fixed=false))
      annotation (Placement(transformation(extent={{24,140},{60,176}})));
    AirMix                         airMix
      annotation (Placement(transformation(extent={{-10,-14},{10,6}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{56,-130},{76,-110}})));
    Modelica.Blocks.Interfaces.RealInput occ "Number of occupants"
      annotation (Placement(transformation(extent={{-180,56},{-152,84}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Interfaces.RealOutput vetot "Total airflow supply [m3]"
      annotation (Placement(transformation(extent={{160,-52},{180,-32}}),
          iconTransformation(extent={{214,-64},{234,-44}})));
    Modelica.Blocks.Continuous.Integrator integrator1(k=1)
      annotation (Placement(transformation(extent={{118,-52},{138,-32}})));
    Modelica.Blocks.Math.Max max1
      annotation (Placement(transformation(extent={{-70,-12},{-50,8}})));
    Modelica.Blocks.Math.Gain ventilation(k=maxVent)
      annotation (Placement(transformation(extent={{-64,-52},{-44,-32}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-64,-110},{-44,-90}})));
    Modelica.Blocks.Interfaces.RealOutput ve "Ventilation airflow rate [m3/h]"
      annotation (Placement(transformation(extent={{160,-88},{180,-68}}),
          iconTransformation(extent={{214,-90},{234,-70}})));
    Modelica.Blocks.Interfaces.RealOutput qr "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{160,-170},{180,-150}}),
          iconTransformation(extent={{214,-160},{234,-140}})));
    Modelica.Blocks.Math.Gain scale(k=0.01) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-104,-42})));
    Modelica.Blocks.Interfaces.RealInput dpos "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-180,-56},{-152,-28}}),
          iconTransformation(extent={{-234,-60},{-206,-32}})));
    Modelica.Blocks.Sources.Constant const2(k=Tve)
      annotation (Placement(transformation(extent={{-152,-18},{-132,2}})));

    Modelica.Blocks.Interfaces.RealInput vpos "Radiator valve position [%]"
      annotation (Placement(transformation(extent={{-180,-114},{-152,-86}}),
          iconTransformation(extent={{-232,-152},{-204,-124}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{16,-14},{36,6}})));
    Modelica.Blocks.Interfaces.RealOutput qv "Ventilation heat load [W]"
      annotation (Placement(transformation(extent={{160,52},{180,72}}),
          iconTransformation(extent={{214,2},{234,22}})));
    Modelica.Blocks.Continuous.Integrator integrator2(k=1)
      annotation (Placement(transformation(extent={{120,12},{140,32}})));
    Modelica.Blocks.Interfaces.RealOutput Qv "Ventilation heat load [Wh]"
      annotation (Placement(transformation(extent={{160,12},{180,32}}),
          iconTransformation(extent={{214,-26},{234,-6}})));
    Modelica.Blocks.Math.Gain scale1(k=0.01)
                                            annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-104,-100})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-166,154},{-138,154}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-115,
            154},{-98,154}},             color={0,0,127}));
    connect(hmltp.y, occHeatGain.Q_flow)
      annotation (Line(points={{-5,64},{6,64}},             color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-14,-100},{4,-100}},                          color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-80,110},{-60,110}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{82,140},{102,140}},            color={0,0,127}));
    connect(T, fromKelvin.Celsius)
      annotation (Line(points={{170,140},{125,140}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-123,110},{-102,110}},         color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-146,110},{-166,110}},         color={0,0,127}));
    connect(hmltp.u2, occeffectiv.y) annotation (Line(points={{-28,58},{-38,58},
            {-38,46},{-43,46}},     color={0,0,127}));
    connect(MetabolicHeat.y[1], occeffectiv.u) annotation (Line(points={{-72.8,
            46},{-66,46}},                color={0,0,127}));
    connect(temperatureSensor.T, MetabolicHeat.u[1]) annotation (Line(points={{82,140},
            {86,140},{86,24},{-110,24},{-110,46},{-100.4,46}},
          color={0,0,127}));
    connect(re.port_b, cair.port) annotation (Line(points={{-40,110},{-20,110},
            {-20,140},{42,140}},color={191,0,0}));
    connect(prescribedHeatFlow.port, cair.port) annotation (Line(points={{-78,154},
            {-20,154},{-20,140},{42,140}},      color={191,0,0}));
    connect(temperatureSensor.port, cair.port)
      annotation (Line(points={{62,140},{42,140}},          color={191,0,0}));
    connect(occHeatGain.port, cair.port)
      annotation (Line(points={{26,64},{42,64},{42,140}},color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{82,140},{
            86,140},{86,-30},{0,-30},{0,-14.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, cair.port) annotation (Line(points={{24,-100},
            {42,-100},{42,140}}, color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{14,-110},
            {14,-120},{54,-120}},                        color={0,0,127}));
    connect(Qr, integrator.y)
      annotation (Line(points={{170,-120},{77,-120}},  color={0,0,127}));
    connect(occ, hmltp.u1) annotation (Line(points={{-166,70},{-28,70}},
                  color={0,0,127}));
    connect(integrator1.y, vetot)
      annotation (Line(points={{139,-42},{170,-42}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, max1.u1) annotation (Line(points={{-123,110},{
            -118,110},{-118,4},{-72,4}},color={0,0,127}));
    connect(max1.y, airMix.Tve) annotation (Line(points={{-49,-2},{-40,-2},{-40,
            -0.6},{-10.8,-0.6}},       color={0,0,127}));
    connect(ventilation.y, airMix.Vve) annotation (Line(points={{-43,-42},{-28,
            -42},{-28,-9},{-11,-9}},   color={0,0,127}));
    connect(ventilation.y, integrator1.u)
      annotation (Line(points={{-43,-42},{116,-42}}, color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-34,
            -100},{-43,-100}},                   color={0,0,127}));
    connect(ventilation.y, ve) annotation (Line(points={{-43,-42},{100.5,-42},{
            100.5,-78},{170,-78}},   color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qr) annotation (Line(points={{14,-110},{14,
            -160},{170,-160}}, color={0,0,127}));
    connect(ventilation.u, scale.y)
      annotation (Line(points={{-66,-42},{-93,-42}},  color={0,0,127}));
    connect(dpos, scale.u)
      annotation (Line(points={{-166,-42},{-116,-42}}, color={0,0,127}));
    connect(max1.u2, toKelvin.Kelvin)
      annotation (Line(points={{-72,-8},{-97,-8}},     color={0,0,127}));
    connect(toKelvin.Celsius, const2.y)
      annotation (Line(points={{-120,-8},{-131,-8}},   color={0,0,127}));
    connect(heatFlowSensor1.port_b, cair.port)
      annotation (Line(points={{36,-4},{42,-4},{42,140}},   color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qv) annotation (Line(points={{26,-14},{26,
            -36},{100,-36},{100,62},{170,62}}, color={0,0,127}));
    connect(integrator2.y, Qv)
      annotation (Line(points={{141,22},{170,22}},   color={0,0,127}));
    connect(integrator2.u, qv) annotation (Line(points={{118,22},{100,22},{100,
            62},{170,62}},     color={0,0,127}));
    connect(vpos, scale1.u)
      annotation (Line(points={{-166,-100},{-116,-100}}, color={0,0,127}));
    connect(scale1.y, heating.u)
      annotation (Line(points={{-93,-100},{-66,-100}},  color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{10,-4},{16,-4}}, color={191,0,0}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-160,
              -180},{160,180}},
          initialScale=0.1)),
      experiment(
        Tolerance=1e-09,
        __Dymola_fixedstepsize=1,
        __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-160,-180},{160,
              180}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end ZoneR1C1_compact;

  model ZoneR1C1vpos "Simple zone"

    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real maxHeat=5000 "Heating power of radiators [W]";
    parameter Real Tve=21 "Ventilation air temperature";
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput T "[K]"
      annotation (Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-120,80},{-100,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{62,110},{82,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Blocks.Math.Product hmltp
      annotation (Placement(transformation(extent={{-62,18},{-42,38}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-14,18},{6,38}})));
    Modelica.Blocks.Tables.CombiTable1D MetabolicHeat(table=[293.15,84.; 325.15,0.])
      annotation (Placement(transformation(extent={{-156,-2},{-132,22}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-146},{-4,-126}})));
    Modelica.Blocks.Interfaces.RealOutput Qr "Heat supplied by radiators [Wh]"
      annotation (Placement(transformation(extent={{216,-170},{236,-150}}),
          iconTransformation(extent={{214,-132},{234,-112}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
      annotation (Placement(transformation(extent={{-166,-54},{-146,-34}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-204,80},{-184,100}})));
    Modelica.Blocks.Math.Gain occeffectiv(k=occheff)
      annotation (Placement(transformation(extent={{-116,0},{-96,20}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cair(
                                    C=tmass*1.2*1005*Vi, T(fixed=false))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    AirMix                         airMix
      annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-170},{170,-150}})));
    Modelica.Blocks.Interfaces.RealInput occ "Number of occupants"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Interfaces.RealOutput vetot "Total airflow supply [m3]"
      annotation (Placement(transformation(extent={{216,-88},{236,-68}}),
          iconTransformation(extent={{214,-64},{234,-44}})));
    Modelica.Blocks.Continuous.Integrator integrator1(k=1)
      annotation (Placement(transformation(extent={{150,-88},{170,-68}})));
    Modelica.Blocks.Math.Max max1
      annotation (Placement(transformation(extent={{-130,-48},{-110,-28}})));
    Modelica.Blocks.Math.Gain ventilation(k=maxVent)
      annotation (Placement(transformation(extent={{-88,-88},{-68,-68}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-88,-146},{-68,-126}})));
    Modelica.Blocks.Interfaces.RealOutput ve "Ventilation airflow rate [m3/h]"
      annotation (Placement(transformation(extent={{216,-128},{236,-108}}),
          iconTransformation(extent={{214,-90},{234,-70}})));
    Modelica.Blocks.Interfaces.RealOutput qr "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-210},{236,-190}}),
          iconTransformation(extent={{214,-160},{234,-140}})));
    Modelica.Blocks.Math.Gain scale(k=0.01) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-150,-78})));
    Modelica.Blocks.Sources.Constant const2(k=Tve)
      annotation (Placement(transformation(extent={{-202,-54},{-182,-34}})));

    Modelica.Blocks.Interfaces.RealInput vpos "Radiator valve position [%]"
      annotation (Placement(transformation(extent={{-240,-150},{-212,-122}}),
          iconTransformation(extent={{-232,-152},{-204,-124}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{16,-60},{36,-40}})));
    Modelica.Blocks.Interfaces.RealOutput qv "Ventilation heat load [W]"
      annotation (Placement(transformation(extent={{216,24},{236,44}}),
          iconTransformation(extent={{214,2},{234,22}})));
    Modelica.Blocks.Continuous.Integrator integrator2(k=1)
      annotation (Placement(transformation(extent={{152,-24},{172,-4}})));
    Modelica.Blocks.Interfaces.RealOutput Qv "Ventilation heat load [Wh]"
      annotation (Placement(transformation(extent={{216,-24},{236,-4}}),
          iconTransformation(extent={{214,-26},{234,-6}})));
    Modelica.Blocks.Math.Gain scale1(k=0.01)
                                            annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-150,-136})));
    Modelica.Blocks.Sources.Constant const1(k=0)
      annotation (Placement(transformation(extent={{-202,-88},{-182,-68}})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(hmltp.y, occHeatGain.Q_flow)
      annotation (Line(points={{-41,28},{-14,28}},          color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-24,-136}},                        color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-120,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{82,120},{180,120}},            color={0,0,127}));
    connect(T, fromKelvin.Celsius)
      annotation (Line(points={{224,120},{210,120},{203,120}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-183,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-206,90},{-226,90}},           color={0,0,127}));
    connect(hmltp.u2, occeffectiv.y) annotation (Line(points={{-64,22},{-80,22},{-80,
            10},{-95,10}},          color={0,0,127}));
    connect(MetabolicHeat.y[1], occeffectiv.u) annotation (Line(points={{-130.8,10},
            {-118,10}},                   color={0,0,127}));
    connect(temperatureSensor.T, MetabolicHeat.u[1]) annotation (Line(points={{82,120},
            {90,120},{90,-14},{-168,-14},{-168,10},{-158.4,10}},
          color={0,0,127}));
    connect(re.port_b, cair.port) annotation (Line(points={{-100,90},{-80,90},{
            -80,120},{46,120}}, color={191,0,0}));
    connect(prescribedHeatFlow.port, cair.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.port, cair.port)
      annotation (Line(points={{62,120},{62,120},{46,120}}, color={191,0,0}));
    connect(occHeatGain.port, cair.port)
      annotation (Line(points={{6,28},{46,28},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{82,120},{
            90,120},{90,-66},{0,-66},{0,-60.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, cair.port) annotation (Line(points={{-4,-136},
            {46,-136},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,-146},
            {-14,-160},{148,-160}},                      color={0,0,127}));
    connect(Qr, integrator.y)
      annotation (Line(points={{226,-160},{171,-160}}, color={0,0,127}));
    connect(occ, hmltp.u1) annotation (Line(points={{-226,34},{-64,34}},
                  color={0,0,127}));
    connect(integrator1.y, vetot)
      annotation (Line(points={{171,-78},{176,-78},{226,-78}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, max1.u1) annotation (Line(points={{-183,90},{-174,
            90},{-174,-32},{-132,-32}}, color={0,0,127}));
    connect(max1.y, airMix.Tve) annotation (Line(points={{-109,-38},{-40,-38},{
            -40,-46.6},{-10.8,-46.6}}, color={0,0,127}));
    connect(ventilation.y, airMix.Vve) annotation (Line(points={{-67,-78},{-28,-78},
            {-28,-55},{-11,-55}},      color={0,0,127}));
    connect(ventilation.y, integrator1.u)
      annotation (Line(points={{-67,-78},{-67,-78},{148,-78}},
                                                     color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,-136},
            {-67,-136}},                         color={0,0,127}));
    connect(ventilation.y, ve) annotation (Line(points={{-67,-78},{126.5,-78},{
            126.5,-118},{226,-118}}, color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qr) annotation (Line(points={{-14,-146},{-14,
            -200},{226,-200}}, color={0,0,127}));
    connect(ventilation.u, scale.y)
      annotation (Line(points={{-90,-78},{-139,-78}}, color={0,0,127}));
    connect(max1.u2, toKelvin.Kelvin)
      annotation (Line(points={{-132,-44},{-145,-44}}, color={0,0,127}));
    connect(toKelvin.Celsius, const2.y)
      annotation (Line(points={{-168,-44},{-181,-44}}, color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{10,-50},{16,-50}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, cair.port)
      annotation (Line(points={{36,-50},{46,-50},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qv) annotation (Line(points={{26,-60},{26,
            -72},{120,-72},{120,34},{226,34}}, color={0,0,127}));
    connect(integrator2.y, Qv)
      annotation (Line(points={{173,-14},{226,-14}}, color={0,0,127}));
    connect(integrator2.u, qv) annotation (Line(points={{150,-14},{120,-14},{
            120,34},{226,34}}, color={0,0,127}));
    connect(vpos, scale1.u)
      annotation (Line(points={{-226,-136},{-162,-136}}, color={0,0,127}));
    connect(scale1.y, heating.u)
      annotation (Line(points={{-139,-136},{-90,-136}}, color={0,0,127}));
    connect(scale.u, const1.y)
      annotation (Line(points={{-162,-78},{-181,-78}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1)),
      experiment(
        Tolerance=1e-09,
        __Dymola_fixedstepsize=1,
        __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end ZoneR1C1vpos;

  model ZoneR2C2 "Simple zone"

    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=1.0 "Internal wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real maxHeat=5000 "Heating power of radiators [W]";
    parameter Real Tve=21 "Ventilation air temperature";
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput T "[K]"
      annotation (Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-120,80},{-100,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{62,110},{82,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Blocks.Math.Product hmltp
      annotation (Placement(transformation(extent={{-62,18},{-42,38}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-14,18},{6,38}})));
    Modelica.Blocks.Tables.CombiTable1D MetabolicHeat(table=[293.15,84.; 325.15,0.])
      annotation (Placement(transformation(extent={{-156,-2},{-132,22}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-146},{-4,-126}})));
    Modelica.Blocks.Interfaces.RealOutput Qr "Heat supplied by radiators [Wh]"
      annotation (Placement(transformation(extent={{216,-170},{236,-150}}),
          iconTransformation(extent={{214,-132},{234,-112}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
      annotation (Placement(transformation(extent={{-166,-54},{-146,-34}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-204,80},{-184,100}})));
    Modelica.Blocks.Math.Gain occeffectiv(k=occheff)
      annotation (Placement(transformation(extent={{-116,0},{-96,20}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cair(
                                    C=tmass*1.2*1005*Vi, T(fixed=false))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    AirMix                         airMix
      annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-170},{170,-150}})));
    Modelica.Blocks.Interfaces.RealInput occ "Number of occupants"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Interfaces.RealOutput vetot "Total airflow supply [m3]"
      annotation (Placement(transformation(extent={{216,-88},{236,-68}}),
          iconTransformation(extent={{214,-64},{234,-44}})));
    Modelica.Blocks.Continuous.Integrator integrator1(k=1)
      annotation (Placement(transformation(extent={{150,-88},{170,-68}})));
    Modelica.Blocks.Math.Max max1
      annotation (Placement(transformation(extent={{-130,-48},{-110,-28}})));
    Modelica.Blocks.Math.Gain ventilation(k=maxVent)
      annotation (Placement(transformation(extent={{-88,-88},{-68,-68}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-88,-146},{-68,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cint(C=imass*1.2*
          1005*Vi, T(fixed=false))
      annotation (Placement(transformation(extent={{108,158},{128,178}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput ve "Ventilation airflow rate [m3/h]"
      annotation (Placement(transformation(extent={{216,-128},{236,-108}}),
          iconTransformation(extent={{214,-90},{234,-70}})));
    Modelica.Blocks.Interfaces.RealOutput qr "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-210},{236,-190}}),
          iconTransformation(extent={{214,-160},{234,-140}})));
    Modelica.Blocks.Math.Gain scale(k=0.01) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-150,-78})));
    Modelica.Blocks.Interfaces.RealInput dpos "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-242,-92},{-214,-64}}),
          iconTransformation(extent={{-234,-60},{-206,-32}})));
    Modelica.Blocks.Sources.Constant const2(k=Tve)
      annotation (Placement(transformation(extent={{-202,-54},{-182,-34}})));

    Modelica.Blocks.Interfaces.RealInput vpos "Radiator valve position [%]"
      annotation (Placement(transformation(extent={{-240,-150},{-212,-122}}),
          iconTransformation(extent={{-232,-152},{-204,-124}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{16,-60},{36,-40}})));
    Modelica.Blocks.Interfaces.RealOutput qv "Ventilation heat load [W]"
      annotation (Placement(transformation(extent={{216,24},{236,44}}),
          iconTransformation(extent={{214,2},{234,22}})));
    Modelica.Blocks.Continuous.Integrator integrator2(k=1)
      annotation (Placement(transformation(extent={{152,-24},{172,-4}})));
    Modelica.Blocks.Interfaces.RealOutput Qv "Ventilation heat load [Wh]"
      annotation (Placement(transformation(extent={{216,-24},{236,-4}}),
          iconTransformation(extent={{214,-26},{234,-6}})));
    Modelica.Blocks.Math.Gain scale1(k=0.01)
                                            annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-150,-136})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(hmltp.y, occHeatGain.Q_flow)
      annotation (Line(points={{-41,28},{-14,28}},          color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-24,-136}},                        color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-120,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{82,120},{180,120}},            color={0,0,127}));
    connect(T, fromKelvin.Celsius)
      annotation (Line(points={{224,120},{210,120},{203,120}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-183,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-206,90},{-226,90}},           color={0,0,127}));
    connect(hmltp.u2, occeffectiv.y) annotation (Line(points={{-64,22},{-84,22},{-84,
            10},{-95,10}},          color={0,0,127}));
    connect(MetabolicHeat.y[1], occeffectiv.u) annotation (Line(points={{-130.8,10},
            {-118,10}},                   color={0,0,127}));
    connect(temperatureSensor.T, MetabolicHeat.u[1]) annotation (Line(points={{82,120},
            {90,120},{90,-14},{-168,-14},{-168,10},{-158.4,10}},
          color={0,0,127}));
    connect(re.port_b, cair.port) annotation (Line(points={{-100,90},{-80,90},{
            -80,120},{46,120}}, color={191,0,0}));
    connect(prescribedHeatFlow.port, cair.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.port, cair.port)
      annotation (Line(points={{62,120},{62,120},{46,120}}, color={191,0,0}));
    connect(occHeatGain.port, cair.port)
      annotation (Line(points={{6,28},{46,28},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{82,120},{
            90,120},{90,-66},{0,-66},{0,-60.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, cair.port) annotation (Line(points={{-4,-136},
            {46,-136},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,-146},
            {-14,-160},{148,-160}},                      color={0,0,127}));
    connect(Qr, integrator.y)
      annotation (Line(points={{226,-160},{171,-160}}, color={0,0,127}));
    connect(occ, hmltp.u1) annotation (Line(points={{-226,34},{-64,34}},
                  color={0,0,127}));
    connect(integrator1.y, vetot)
      annotation (Line(points={{171,-78},{176,-78},{226,-78}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, max1.u1) annotation (Line(points={{-183,90},{-174,
            90},{-174,-32},{-132,-32}}, color={0,0,127}));
    connect(max1.y, airMix.Tve) annotation (Line(points={{-109,-38},{-40,-38},{
            -40,-46.6},{-10.8,-46.6}}, color={0,0,127}));
    connect(ventilation.y, airMix.Vve) annotation (Line(points={{-67,-78},{-28,-78},
            {-28,-55},{-11,-55}},      color={0,0,127}));
    connect(ventilation.y, integrator1.u)
      annotation (Line(points={{-67,-78},{-67,-78},{148,-78}},
                                                     color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,-136},
            {-67,-136}},                         color={0,0,127}));
    connect(re1.port_b, cint.port) annotation (Line(points={{96,158},{110,158},
            {118,158}}, color={191,0,0}));
    connect(re1.port_a, cair.port) annotation (Line(points={{76,158},{62,158},{
            62,120},{46,120}}, color={191,0,0}));
    connect(ventilation.y, ve) annotation (Line(points={{-67,-78},{126.5,-78},{
            126.5,-118},{226,-118}}, color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qr) annotation (Line(points={{-14,-146},{-14,
            -200},{226,-200}}, color={0,0,127}));
    connect(ventilation.u, scale.y)
      annotation (Line(points={{-90,-78},{-139,-78}}, color={0,0,127}));
    connect(dpos, scale.u)
      annotation (Line(points={{-228,-78},{-162,-78}}, color={0,0,127}));
    connect(max1.u2, toKelvin.Kelvin)
      annotation (Line(points={{-132,-44},{-145,-44}}, color={0,0,127}));
    connect(toKelvin.Celsius, const2.y)
      annotation (Line(points={{-168,-44},{-181,-44}}, color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{10,-50},{16,-50}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, cair.port)
      annotation (Line(points={{36,-50},{46,-50},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qv) annotation (Line(points={{26,-60},{26,
            -72},{120,-72},{120,34},{226,34}}, color={0,0,127}));
    connect(integrator2.y, Qv)
      annotation (Line(points={{173,-14},{226,-14}}, color={0,0,127}));
    connect(integrator2.u, qv) annotation (Line(points={{150,-14},{120,-14},{
            120,34},{226,34}}, color={0,0,127}));
    connect(vpos, scale1.u)
      annotation (Line(points={{-226,-136},{-162,-136}}, color={0,0,127}));
    connect(scale1.y, heating.u)
      annotation (Line(points={{-139,-136},{-90,-136}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1)),
      experiment(
        Tolerance=1e-09,
        __Dymola_fixedstepsize=1,
        __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end ZoneR2C2;

  model ZoneCO2R1C1 "Simple zone"
    import MShootBS2019;

    // Temperature related parameters
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real maxHeat=5000 "Heating power of radiators [W]";
    parameter Real Tve=21 "Ventilation air temperature";
    parameter Real occheff=1. "Occupant heat generation effectiveness";

    // CO2 related parameters
    parameter Real CO2pp=0.02 "CO2 generation per person [m3/h]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real Vinf=100 "Infiltration flow rate [m3/h]";

    // Parameters affecting both CO2 and temperatures
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput T "[K]"
      annotation (Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-120,80},{-100,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{62,110},{82,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Blocks.Math.Product hmltp
      annotation (Placement(transformation(extent={{-62,18},{-42,38}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-14,18},{6,38}})));
    Modelica.Blocks.Tables.CombiTable1D MetabolicHeat(table=[293.15,84.; 325.15,0.])
      annotation (Placement(transformation(extent={{-156,-2},{-132,22}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-146},{-4,-126}})));
    Modelica.Blocks.Interfaces.RealOutput Qr "Heat supplied by radiators [Wh]"
      annotation (Placement(transformation(extent={{216,-170},{236,-150}}),
          iconTransformation(extent={{214,-132},{234,-112}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
      annotation (Placement(transformation(extent={{-166,-54},{-146,-34}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-204,80},{-184,100}})));
    Modelica.Blocks.Math.Gain occeffectiv(k=occheff)
      annotation (Placement(transformation(extent={{-116,0},{-96,20}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cair(
                                    C=tmass*1.2*1005*Vi, T(fixed=false))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    AirMix                         airMix
      annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-170},{170,-150}})));
    Modelica.Blocks.Interfaces.RealInput occ "Number of occupants"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Interfaces.RealOutput vetot "Total airflow supply [m3]"
      annotation (Placement(transformation(extent={{216,-88},{236,-68}}),
          iconTransformation(extent={{214,-64},{234,-44}})));
    Modelica.Blocks.Continuous.Integrator integrator1(k=1)
      annotation (Placement(transformation(extent={{150,-88},{170,-68}})));
    Modelica.Blocks.Math.Max max1
      annotation (Placement(transformation(extent={{-130,-48},{-110,-28}})));
    Modelica.Blocks.Math.Gain ventilation(k=maxVent)
      annotation (Placement(transformation(extent={{-88,-88},{-68,-68}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-88,-146},{-68,-126}})));
    Modelica.Blocks.Interfaces.RealOutput ve "Ventilation airflow rate [m3/h]"
      annotation (Placement(transformation(extent={{216,-128},{236,-108}}),
          iconTransformation(extent={{214,-90},{234,-70}})));
    Modelica.Blocks.Interfaces.RealOutput qr "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-210},{236,-190}}),
          iconTransformation(extent={{214,-160},{234,-140}})));
    Modelica.Blocks.Math.Gain scale(k=0.01) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-150,-78})));
    Modelica.Blocks.Interfaces.RealInput dpos "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-242,-92},{-214,-64}}),
          iconTransformation(extent={{-234,-60},{-206,-32}})));
    Modelica.Blocks.Sources.Constant const2(k=Tve)
      annotation (Placement(transformation(extent={{-202,-54},{-182,-34}})));

    Modelica.Blocks.Interfaces.RealInput vpos "Radiator valve position [%]"
      annotation (Placement(transformation(extent={{-240,-150},{-212,-122}}),
          iconTransformation(extent={{-232,-152},{-204,-124}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{16,-60},{36,-40}})));
    Modelica.Blocks.Interfaces.RealOutput qv "Ventilation heat load [W]"
      annotation (Placement(transformation(extent={{216,24},{236,44}}),
          iconTransformation(extent={{214,2},{234,22}})));
    Modelica.Blocks.Continuous.Integrator integrator2(k=1)
      annotation (Placement(transformation(extent={{152,-24},{172,-4}})));
    Modelica.Blocks.Interfaces.RealOutput Qv "Ventilation heat load [Wh]"
      annotation (Placement(transformation(extent={{216,-24},{236,-4}}),
          iconTransformation(extent={{214,-26},{234,-6}})));
    Modelica.Blocks.Math.Gain scale1(k=0.01)
                                            annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-150,-136})));
    MShootBS2019.CO2 co2(
      CO2PerPerson=CO2pp,
      Vi=Vi,
      CO2Neutral=CO2n)
      annotation (Placement(transformation(extent={{16,60},{36,80}})));
    Modelica.Blocks.Sources.Constant const1(k=Vinf)
      annotation (Placement(transformation(extent={{-54,86},{-34,106}})));

    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{-20,70},{0,90}})));
    Modelica.Blocks.Interfaces.RealOutput CO2 "CO2 concentration [ppm]"
      annotation (Placement(transformation(extent={{216,60},{236,80}}),
          iconTransformation(extent={{214,2},{234,22}})));

  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(hmltp.y, occHeatGain.Q_flow)
      annotation (Line(points={{-41,28},{-14,28}},          color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-24,-136}},                        color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-120,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{82,120},{180,120}},            color={0,0,127}));
    connect(T, fromKelvin.Celsius)
      annotation (Line(points={{224,120},{210,120},{203,120}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-183,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-206,90},{-226,90}},           color={0,0,127}));
    connect(hmltp.u2, occeffectiv.y) annotation (Line(points={{-64,22},{-80,22},{-80,
            10},{-95,10}},          color={0,0,127}));
    connect(MetabolicHeat.y[1], occeffectiv.u) annotation (Line(points={{-130.8,10},
            {-118,10}},                   color={0,0,127}));
    connect(temperatureSensor.T, MetabolicHeat.u[1]) annotation (Line(points={{82,120},
            {90,120},{90,-14},{-168,-14},{-168,10},{-158.4,10}},
          color={0,0,127}));
    connect(re.port_b, cair.port) annotation (Line(points={{-100,90},{-80,90},{
            -80,120},{46,120}}, color={191,0,0}));
    connect(prescribedHeatFlow.port, cair.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.port, cair.port)
      annotation (Line(points={{62,120},{62,120},{46,120}}, color={191,0,0}));
    connect(occHeatGain.port, cair.port)
      annotation (Line(points={{6,28},{46,28},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{82,120},{
            90,120},{90,-66},{0,-66},{0,-60.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, cair.port) annotation (Line(points={{-4,-136},
            {46,-136},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,-146},
            {-14,-160},{148,-160}},                      color={0,0,127}));
    connect(Qr, integrator.y)
      annotation (Line(points={{226,-160},{171,-160}}, color={0,0,127}));
    connect(occ, hmltp.u1) annotation (Line(points={{-226,34},{-64,34}},
                  color={0,0,127}));
    connect(integrator1.y, vetot)
      annotation (Line(points={{171,-78},{176,-78},{226,-78}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, max1.u1) annotation (Line(points={{-183,90},{-174,
            90},{-174,-32},{-132,-32}}, color={0,0,127}));
    connect(max1.y, airMix.Tve) annotation (Line(points={{-109,-38},{-40,-38},{
            -40,-46.6},{-10.8,-46.6}}, color={0,0,127}));
    connect(ventilation.y, airMix.Vve) annotation (Line(points={{-67,-78},{-28,-78},
            {-28,-55},{-11,-55}},      color={0,0,127}));
    connect(ventilation.y, integrator1.u)
      annotation (Line(points={{-67,-78},{-67,-78},{148,-78}},
                                                     color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,-136},
            {-67,-136}},                         color={0,0,127}));
    connect(ventilation.y, ve) annotation (Line(points={{-67,-78},{120.5,-78},{120.5,
            -118},{226,-118}},       color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qr) annotation (Line(points={{-14,-146},{-14,
            -200},{226,-200}}, color={0,0,127}));
    connect(ventilation.u, scale.y)
      annotation (Line(points={{-90,-78},{-139,-78}}, color={0,0,127}));
    connect(dpos, scale.u)
      annotation (Line(points={{-228,-78},{-162,-78}}, color={0,0,127}));
    connect(max1.u2, toKelvin.Kelvin)
      annotation (Line(points={{-132,-44},{-145,-44}}, color={0,0,127}));
    connect(toKelvin.Celsius, const2.y)
      annotation (Line(points={{-168,-44},{-181,-44}}, color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{10,-50},{16,-50}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, cair.port)
      annotation (Line(points={{36,-50},{46,-50},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qv) annotation (Line(points={{26,-60},{26,
            -72},{120,-72},{120,34},{226,34}}, color={0,0,127}));
    connect(integrator2.y, Qv)
      annotation (Line(points={{173,-14},{226,-14}}, color={0,0,127}));
    connect(integrator2.u, qv) annotation (Line(points={{150,-14},{120,-14},{
            120,34},{226,34}}, color={0,0,127}));
    connect(vpos, scale1.u)
      annotation (Line(points={{-226,-136},{-162,-136}}, color={0,0,127}));
    connect(scale1.y, heating.u)
      annotation (Line(points={{-139,-136},{-90,-136}}, color={0,0,127}));
    connect(occ, co2.persons) annotation (Line(points={{-226,34},{-80,34},{-80,64},
            {15.6,64}}, color={0,0,127}));
    connect(add.y, co2.Vve) annotation (Line(points={{1,80},{7.5,80},{7.5,76},{15.6,
            76}}, color={0,0,127}));
    connect(const1.y, add.u1) annotation (Line(points={{-33,96},{-28,96},{-28,86},
            {-22,86}}, color={0,0,127}));
    connect(ventilation.y, add.u2) annotation (Line(points={{-67,-78},{-28,-78},{-28,
            74},{-22,74}}, color={0,0,127}));
    connect(co2.CO2, CO2)
      annotation (Line(points={{36.6,70},{226,70}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1)),
      experiment(
        Tolerance=1e-09,
        __Dymola_fixedstepsize=1,
        __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end ZoneCO2R1C1;

  model AirMix
    "Calculates heat flow based on indoor temperature and ventilation airflow and temperature"

    Modelica.Blocks.Sources.Constant h_to_s(k=1/3600)
      annotation (Placement(transformation(extent={{-64,58},{-44,78}})));
    Modelica.Blocks.Sources.Constant AirRhoCp(k=1.2*1005)
      annotation (Placement(transformation(extent={{-64,26},{-44,46}})));
    Modelica.Blocks.Math.Add add(k1=-1, k2=1) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-42,4})));
    Modelica.Blocks.Math.MultiProduct multiProduct(      significantDigits=8, nu=4)
      annotation (Placement(transformation(extent={{-10,22},{14,46}})));
    Modelica.Blocks.Interfaces.RealInput Tve
      annotation (Placement(transformation(extent={{-128,14},{-88,54}})));
    Modelica.Blocks.Interfaces.RealInput Ti annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=90,
          origin={0,-108})));
    Modelica.Blocks.Interfaces.RealInput Vve "[m3/h]"
      annotation (Placement(transformation(extent={{-130,-70},{-90,-30}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{48,-10},{68,10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  equation
    connect(h_to_s.y, multiProduct.u[1]) annotation (Line(points={{-43,68},{-22,
            68},{-22,40.3},{-10,40.3}},
                                      color={0,0,127}));
    connect(AirRhoCp.y, multiProduct.u[2]) annotation (Line(points={{-43,36},{
            -43,36.1},{-10,36.1}},
                                 color={0,0,127}));
    connect(add.y, multiProduct.u[3]) annotation (Line(points={{-31,4},{-22,4},
            {-22,31.9},{-10,31.9}},
                           color={0,0,127}));
    connect(Tve, add.u2) annotation (Line(points={{-108,34},{-80,34},{-80,10},{-54,
            10}}, color={0,0,127}));
    connect(Ti, add.u1) annotation (Line(points={{0,-108},{0,-26},{-68,-26},{-68,-2},
            {-54,-2}}, color={0,0,127}));
    connect(Vve, multiProduct.u[4]) annotation (Line(points={{-110,-50},{-18,
            -50},{-18,27.7},{-10,27.7}},
                                      color={0,0,127}));
    connect(multiProduct.y, prescribedHeatFlow.Q_flow) annotation (Line(points=
            {{16.04,34},{32,34},{32,0},{48,0}}, color={0,0,127}));
    connect(prescribedHeatFlow.port, port_b)
      annotation (Line(points={{68,0},{74,0},{100,0}}, color={191,0,0}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
            extent={{-100,-100},{100,100}}), graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={28,108,200},
            lineThickness=0.5),
          Line(points={{-68,46},{-62,28},{-56,46}}, color={28,108,200}),
          Line(points={{-78,46},{-78,26}}, color={28,108,200}),
          Line(points={{-74,46},{-82,46}}, color={28,108,200}),
          Line(points={{-84,-38},{-78,-56},{-72,-38}}, color={28,108,200}),
          Line(points={{-68,-38},{-62,-56},{-56,-38}}, color={28,108,200}),
          Line(points={{-8,-62},{-8,-82}}, color={28,108,200}),
          Line(points={{-4,-62},{-12,-62}}, color={28,108,200}),
          Line(points={{2,-62},{2,-82}}, color={28,108,200}),
          Line(points={{0,-50},{0,0},{90,0}}, color={28,108,200}),
          Line(points={{-50,36},{-40,36},{-20,0},{0,0}}, color={28,108,200}),
          Line(points={{-50,-48},{-40,-48},{-20,0}}, color={28,108,200}),
          Ellipse(
            extent={{-4,4},{4,-4}},
            lineColor={28,108,200},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid)}));
  end AirMix;

  model CO2Eqn "CO2 balance model (equations)"

    Modelica.Blocks.Interfaces.RealInput Vve
      "Ventilation air flow rate [m3/h]"
      annotation (Placement(transformation(extent={{-116,62},{-92,86}})));
    Modelica.Blocks.Interfaces.RealInput CO2ppmv_s
      "CO2 concentration in incoming air [ppmv]"
      annotation (Placement(transformation(extent={{-116,14},{-92,38}})));
    Modelica.Blocks.Interfaces.RealInput persons "Number of persons"
      annotation (Placement(transformation(extent={{-116,-38},{-92,-14}})));

    Modelica.Blocks.Interfaces.RealInput CO2_per_person
      "CO2 generation per person [ppmv/h]"
      annotation (Placement(transformation(extent={{-116,-86},{-92,-62}})));
    Modelica.Blocks.Interfaces.RealOutput CO2output
      annotation (Placement(transformation(extent={{96,-10},{116,10}})));

    parameter Real Vi "Indoor volume [m3]";

    // Real VCO2_g_rate "CO2 human generation rate [m3/h]";
    // Real VCO2_s_rate "CO2 supply rate [m3/h]";
    // Real VCO2_ex_rate "C02 extraction rate [m3/h]";
    // Real VCO2_i "CO2 indoor concentration [m3]";
    Real CO2ppmv_i "CO2 indoor concentration [ppmv]";

  equation
    // Short version of the balance (works with EstimationPy!):
    3600 * Vi * der(CO2ppmv_i) = persons * CO2_per_person * 10^6 + Vve * (CO2ppmv_s - CO2ppmv_i);
    CO2output = CO2ppmv_i;

    // Long version (more intuitive) of the balance:
    /*
  der(VCO2_i) = VCO2_g_rate + VCO2_s_rate - VCO2_ex_rate "Transient balance";
  VCO2_g_rate = persons * CO2_per_person / 3600 "Generation rate";
  VCO2_s_rate = Vve * CO2ppmv_s / 10^6 / 3600 "Supply rate";
  VCO2_ex_rate = Vve * CO2ppmv_i / 10^6 / 3600 "Extraction rate";
  VCO2_i = CO2ppmv_i * Vi / 10^6 "Relationship between PPMV and m3";
  CO2output = CO2ppmv_i "Output in m3";
  */

     annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
            extent={{-100,-100},{100,100}}), graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Line(
            points={{84,-4},{-16,-44}},
            color={28,108,200},
            thickness=0.5),
          Polygon(
            points={{34,-24},{26,-44},{42,-44},{34,-24}},
            lineColor={28,108,200},
            lineThickness=0.5,
            fillPattern=FillPattern.Sphere,
            fillColor={28,108,200}),
          Line(points={{-76,-18},{-56,-18}}, color={255,0,0}),
          Line(points={{-60,-48},{-66,-30},{-66,-10}}, color={255,0,0}),
          Ellipse(extent={{-72,2},{-60,-10}}, lineColor={255,0,0}),
          Line(points={{-66,-30},{-74,-48}}, color={255,0,0}),
          Ellipse(
            extent={{-76,-84},{-74,-86}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-72,-80},{-70,-82}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-64,-78},{-62,-80}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-68,-72},{-66,-74}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-60,-62},{-58,-64}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-74,-64},{-72,-66}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-76,-74},{-74,-76}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-88,24},{-82,8},{-72,28},{-62,8},{-52,28},{-44,14}},
            color={0,255,255},
            thickness=0.5,
            smooth=Smooth.Bezier),
          Line(
            points={{-84,40},{-78,24},{-68,44},{-58,24},{-48,44},{-40,30}},
            color={0,255,255},
            thickness=0.5,
            smooth=Smooth.Bezier),
          Ellipse(
            extent={{-54,18},{-52,16}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-56,36},{-54,34}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-78,36},{-76,34}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-72,30},{-70,28}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-88,70},{-82,54},{-72,74},{-62,54},{-52,74},{-44,60}},
            color={0,255,255},
            thickness=0.5,
            smooth=Smooth.Bezier),
          Line(
            points={{-86,78},{-80,62},{-70,82},{-60,62},{-50,82},{-42,68}},
            color={0,255,255},
            thickness=0.5,
            smooth=Smooth.Bezier),
          Line(
            points={{-84,86},{-78,70},{-68,90},{-58,70},{-48,90},{-40,76}},
            color={0,255,255},
            thickness=0.5,
            smooth=Smooth.Bezier),
          Ellipse(
            extent={{-64,26},{-62,24}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-74,16},{-72,14}},
            lineColor={135,135,135},
            lineThickness=0.5,
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Line(
            points={{24,38},{10,38},{8,36},{8,0},{10,-2},{24,-2}},
            color={0,0,127},
            thickness=0.5),
          Line(
            points={{46,38},{30,38},{28,36},{28,0},{30,-2},{46,-2},{48,0},{
                48,36},{46,38}},
            color={0,0,127},
            thickness=0.5),
          Line(
            points={{52,2},{54,4},{56,4},{58,2},{58,0},{52,-6},{58,-6}},
            color={0,0,127},
            thickness=0.5)}));
  end CO2Eqn;

  model CO2 "CO2 balance model"

    Modelica.Blocks.Interfaces.RealInput Vve "[m3/h]"
      annotation (Placement(transformation(extent={{-120,44},{-88,76}}),
          iconTransformation(extent={{-120,44},{-88,76}})));
    Modelica.Blocks.Interfaces.RealOutput CO2
      annotation (Placement(transformation(extent={{96,-10},{116,10}})));
    CO2Eqn balance(Vi=Vi)
      annotation (Placement(transformation(extent={{-22,-25},{26,25}})));
    Modelica.Blocks.Interfaces.RealInput persons
      annotation (Placement(transformation(extent={{-119,-75},{-89,-45}}),
          iconTransformation(extent={{-119,-75},{-89,-45}})));
    Modelica.Blocks.Sources.Constant CO2Generation(k=CO2PerPerson) "[m3/h]"
      annotation (Placement(transformation(extent={{-76,-38},{-56,-18}})));
    parameter Real CO2PerPerson=0.02 "CO2 generation per person [m3/h]";
    parameter Real Vi=100 "Indoor volume [m3]";
    parameter Real CO2Neutral=400 "CO2 Neutral Level";
    Modelica.Blocks.Sources.Constant CO2NeutralLevel(k=CO2Neutral)
      annotation (Placement(transformation(extent={{-76,6},{-56,26}})));

  equation
    connect(persons, balance.persons) annotation (Line(points={{-104,-60},{
            -84,-60},{-84,-6.5},{-22.96,-6.5}}, color={0,0,127}));
    connect(Vve, balance.Vve) annotation (Line(points={{-104,60},{-32,60},{
            -32,18.5},{-22.96,18.5}}, color={0,0,127}));
    connect(balance.CO2output, CO2) annotation (Line(points={{27.44,0},{
            27.44,0},{106,0}}, color={0,0,127}));
    connect(CO2Generation.y, balance.CO2_per_person) annotation (Line(
          points={{-55,-28},{-40,-28},{-40,-18.5},{-22.96,-18.5}}, color={0,
            0,127}));
    connect(CO2NeutralLevel.y, balance.CO2ppmv_s) annotation (Line(points={
            {-55,16},{-40,16},{-40,6.5},{-22.96,6.5}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
            extent={{-100,-100},{100,100}}), graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200},
            lineThickness=0.5),
          Line(
            points={{36,60}},
            color={28,108,200},
            thickness=0.5),
          Line(
            points={{-66,14}},
            color={135,135,135},
            thickness=0.5),
          Line(
            points={{82,10},{-18,-30}},
            color={28,108,200},
            thickness=0.5),
          Polygon(
            points={{32,-10},{24,-30},{40,-30},{32,-10}},
            lineColor={28,108,200},
            lineThickness=0.5,
            fillPattern=FillPattern.Sphere,
            fillColor={28,108,200}),
          Line(points={{-72,-54},{-52,-54}}, color={255,0,0}),
          Line(points={{-54,-84},{-62,-66},{-62,-46}}, color={255,0,0}),
          Ellipse(extent={{-68,-34},{-56,-46}}, lineColor={255,0,0}),
          Line(points={{-62,-66},{-70,-84}}, color={255,0,0}),
          Line(
            points={{-82,56},{-76,40},{-66,60},{-56,40},{-46,60},{-38,46}},
            color={0,128,255},
            thickness=0.5,
            smooth=Smooth.Bezier),
          Line(
            points={{-78,72},{-72,56},{-62,76},{-52,56},{-42,76},{-34,62}},
            color={0,128,255},
            thickness=0.5,
            smooth=Smooth.Bezier),
          Line(
            points={{22,52},{8,52},{6,50},{6,14},{8,12},{22,12}},
            color={0,0,127},
            thickness=0.5),
          Line(
            points={{44,52},{28,52},{26,50},{26,14},{28,12},{44,12},{46,14},
                {46,50},{44,52}},
            color={0,0,127},
            thickness=0.5),
          Line(
            points={{50,16},{52,18},{54,18},{56,16},{56,14},{50,8},{56,8}},
            color={0,0,127},
            thickness=0.5)}));
  end CO2;

  model ZoneCO2R1C1PID "Simple zone"
    import MShootBS2019;

    // Temperature related parameters
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real maxHeat=5000 "Heating power of radiators [W]";
    parameter Real Tve=21 "Ventilation air temperature";
    parameter Real occheff=1. "Occupant heat generation effectiveness";

    // CO2 related parameters
    parameter Real CO2pp=0.02 "CO2 generation per person [m3/h]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real Vinf=100 "Infiltration flow rate [m3/h]";

    // Parameters affecting both CO2 and temperatures
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";

    // Parameters related to the temp. and co2 controllers
    parameter Real kco2=10 "Gain of CO2 controller";
    parameter Real kt=100 "Gain of temperature controller";

    Modelica.Blocks.Interfaces.RealInput solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput T "[K]"
      annotation (Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-120,80},{-100,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{62,110},{82,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Blocks.Math.Product hmltp
      annotation (Placement(transformation(extent={{-62,18},{-42,38}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-14,18},{6,38}})));
    Modelica.Blocks.Tables.CombiTable1D MetabolicHeat(table=[293.15,84.; 325.15,0.])
      annotation (Placement(transformation(extent={{-156,-2},{-132,22}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-146},{-4,-126}})));
    Modelica.Blocks.Interfaces.RealOutput Qr "Heat supplied by radiators [Wh]"
      annotation (Placement(transformation(extent={{216,-170},{236,-150}}),
          iconTransformation(extent={{214,-132},{234,-112}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin
      annotation (Placement(transformation(extent={{-166,-54},{-146,-34}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-204,80},{-184,100}})));
    Modelica.Blocks.Math.Gain occeffectiv(k=occheff)
      annotation (Placement(transformation(extent={{-116,0},{-96,20}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cair(
                                    C=tmass*1.2*1005*Vi, T(fixed=false))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    AirMix                         airMix
      annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-170},{170,-150}})));
    Modelica.Blocks.Interfaces.RealInput occ "Number of occupants"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Interfaces.RealOutput vetot "Total airflow supply [m3]"
      annotation (Placement(transformation(extent={{216,-88},{236,-68}}),
          iconTransformation(extent={{214,-64},{234,-44}})));
    Modelica.Blocks.Continuous.Integrator integrator1(k=1)
      annotation (Placement(transformation(extent={{150,-88},{170,-68}})));
    Modelica.Blocks.Math.Max max1
      annotation (Placement(transformation(extent={{-130,-48},{-110,-28}})));
    Modelica.Blocks.Math.Gain ventilation(k=maxVent)
      annotation (Placement(transformation(extent={{-82,-88},{-62,-68}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-84,-146},{-64,-126}})));
    Modelica.Blocks.Interfaces.RealOutput ve "Ventilation airflow rate [m3/h]"
      annotation (Placement(transformation(extent={{216,-128},{236,-108}}),
          iconTransformation(extent={{214,-90},{234,-70}})));
    Modelica.Blocks.Interfaces.RealOutput qr "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-210},{236,-190}}),
          iconTransformation(extent={{214,-160},{234,-140}})));
    Modelica.Blocks.Math.Gain scale(k=-0.01)
                                            annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-106,-78})));
    Modelica.Blocks.Sources.Constant const2(k=Tve)
      annotation (Placement(transformation(extent={{-202,-54},{-182,-34}})));

    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{16,-60},{36,-40}})));
    Modelica.Blocks.Interfaces.RealOutput qv "Ventilation heat load [W]"
      annotation (Placement(transformation(extent={{216,24},{236,44}}),
          iconTransformation(extent={{214,2},{234,22}})));
    Modelica.Blocks.Continuous.Integrator integrator2(k=1)
      annotation (Placement(transformation(extent={{152,-24},{172,-4}})));
    Modelica.Blocks.Interfaces.RealOutput Qv "Ventilation heat load [Wh]"
      annotation (Placement(transformation(extent={{216,-24},{236,-4}}),
          iconTransformation(extent={{214,-26},{234,-6}})));
    Modelica.Blocks.Math.Gain scale1(k=0.01)
                                            annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-106,-136})));
    MShootBS2019.CO2 co2(
      CO2PerPerson=CO2pp,
      Vi=Vi,
      CO2Neutral=CO2n)
      annotation (Placement(transformation(extent={{16,60},{36,80}})));
    Modelica.Blocks.Sources.Constant const1(k=Vinf)
      annotation (Placement(transformation(extent={{-54,86},{-34,106}})));

    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{-20,70},{0,90}})));
    Modelica.Blocks.Interfaces.RealOutput CO2 "CO2 concentration [ppm]"
      annotation (Placement(transformation(extent={{216,60},{236,80}}),
          iconTransformation(extent={{214,2},{234,22}})));

    Modelica.Blocks.Continuous.LimPID PIDt(
      yMax=100,
      yMin=0,
      k=kt,
      controllerType=Modelica.Blocks.Types.SimpleController.PI)
      annotation (Placement(transformation(extent={{-166,-146},{-146,-126}})));

    Modelica.Blocks.Interfaces.RealOutput vpos annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={-138,-226}), iconTransformation(extent={{214,-64},{234,-44}})));
    Modelica.Blocks.Interfaces.RealOutput dpos annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={-70,-228}), iconTransformation(extent={{214,-64},{234,-44}})));
    Modelica.Blocks.Math.Gain scale2(k=-1)  annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-100,-192})));
    Modelica.Blocks.Continuous.LimPID PIDco2(
      k=kco2,
      yMax=0,
      yMin=-100,
      controllerType=Modelica.Blocks.Types.SimpleController.PI)
      annotation (Placement(transformation(extent={{-166,-88},{-146,-68}})));

    Modelica.Blocks.Interfaces.RealInput tstp "Temperature setpoint [K]"
      annotation (Placement(transformation(extent={{-240,-150},{-212,-122}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Interfaces.RealInput co2stp "CO2 setpoint" annotation (
        Placement(transformation(extent={{-240,-92},{-212,-64}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
  equation
    connect(solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(hmltp.y, occHeatGain.Q_flow)
      annotation (Line(points={{-41,28},{-14,28}},          color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-24,-136}},                        color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-120,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{82,120},{180,120}},            color={0,0,127}));
    connect(T, fromKelvin.Celsius)
      annotation (Line(points={{224,120},{210,120},{203,120}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-183,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-206,90},{-226,90}},           color={0,0,127}));
    connect(hmltp.u2, occeffectiv.y) annotation (Line(points={{-64,22},{-80,22},{-80,
            10},{-95,10}},          color={0,0,127}));
    connect(MetabolicHeat.y[1], occeffectiv.u) annotation (Line(points={{-130.8,10},
            {-118,10}},                   color={0,0,127}));
    connect(temperatureSensor.T, MetabolicHeat.u[1]) annotation (Line(points={{82,120},
            {90,120},{90,-14},{-168,-14},{-168,10},{-158.4,10}},
          color={0,0,127}));
    connect(re.port_b, cair.port) annotation (Line(points={{-100,90},{-80,90},{
            -80,120},{46,120}}, color={191,0,0}));
    connect(prescribedHeatFlow.port, cair.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.port, cair.port)
      annotation (Line(points={{62,120},{62,120},{46,120}}, color={191,0,0}));
    connect(occHeatGain.port, cair.port)
      annotation (Line(points={{6,28},{46,28},{46,120}}, color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{82,120},{
            90,120},{90,-66},{0,-66},{0,-60.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, cair.port) annotation (Line(points={{-4,-136},
            {46,-136},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,-146},
            {-14,-160},{148,-160}},                      color={0,0,127}));
    connect(Qr, integrator.y)
      annotation (Line(points={{226,-160},{171,-160}}, color={0,0,127}));
    connect(occ, hmltp.u1) annotation (Line(points={{-226,34},{-64,34}},
                  color={0,0,127}));
    connect(integrator1.y, vetot)
      annotation (Line(points={{171,-78},{176,-78},{226,-78}},
                                                     color={0,0,127}));
    connect(toKelvin1.Kelvin, max1.u1) annotation (Line(points={{-183,90},{-174,
            90},{-174,-32},{-132,-32}}, color={0,0,127}));
    connect(max1.y, airMix.Tve) annotation (Line(points={{-109,-38},{-40,-38},{
            -40,-46.6},{-10.8,-46.6}}, color={0,0,127}));
    connect(ventilation.y, airMix.Vve) annotation (Line(points={{-61,-78},{-28,-78},
            {-28,-55},{-11,-55}},      color={0,0,127}));
    connect(ventilation.y, integrator1.u)
      annotation (Line(points={{-61,-78},{148,-78}}, color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,-136},
            {-63,-136}},                         color={0,0,127}));
    connect(ventilation.y, ve) annotation (Line(points={{-61,-78},{120.5,-78},{120.5,
            -118},{226,-118}},       color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qr) annotation (Line(points={{-14,-146},{-14,
            -200},{226,-200}}, color={0,0,127}));
    connect(ventilation.u, scale.y)
      annotation (Line(points={{-84,-78},{-95,-78}},  color={0,0,127}));
    connect(max1.u2, toKelvin.Kelvin)
      annotation (Line(points={{-132,-44},{-145,-44}}, color={0,0,127}));
    connect(toKelvin.Celsius, const2.y)
      annotation (Line(points={{-168,-44},{-181,-44}}, color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{10,-50},{16,-50}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, cair.port)
      annotation (Line(points={{36,-50},{46,-50},{46,120}}, color={191,0,0}));
    connect(heatFlowSensor1.Q_flow, qv) annotation (Line(points={{26,-60},{26,
            -72},{120,-72},{120,34},{226,34}}, color={0,0,127}));
    connect(integrator2.y, Qv)
      annotation (Line(points={{173,-14},{226,-14}}, color={0,0,127}));
    connect(integrator2.u, qv) annotation (Line(points={{150,-14},{120,-14},{
            120,34},{226,34}}, color={0,0,127}));
    connect(scale1.y, heating.u)
      annotation (Line(points={{-95,-136},{-86,-136}},  color={0,0,127}));
    connect(occ, co2.persons) annotation (Line(points={{-226,34},{-80,34},{-80,64},
            {15.6,64}}, color={0,0,127}));
    connect(add.y, co2.Vve) annotation (Line(points={{1,80},{7.5,80},{7.5,76},{15.6,
            76}}, color={0,0,127}));
    connect(const1.y, add.u1) annotation (Line(points={{-33,96},{-28,96},{-28,86},
            {-22,86}}, color={0,0,127}));
    connect(ventilation.y, add.u2) annotation (Line(points={{-61,-78},{-28,-78},{-28,
            74},{-22,74}}, color={0,0,127}));
    connect(co2.CO2, CO2)
      annotation (Line(points={{36.6,70},{226,70}}, color={0,0,127}));
    connect(scale1.u, PIDt.y)
      annotation (Line(points={{-118,-136},{-145,-136}}, color={0,0,127}));
    connect(temperatureSensor.T, PIDt.u_m) annotation (Line(points={{82,120},{90,120},
            {90,-172},{-156,-172},{-156,-148}}, color={0,0,127}));
    connect(PIDt.y, vpos) annotation (Line(points={{-145,-136},{-138,-136},{-138,-226}},
          color={0,0,127}));
    connect(scale2.y, dpos) annotation (Line(points={{-89,-192},{-70,-192},{-70,-228}},
          color={0,0,127}));
    connect(scale.u, PIDco2.y)
      annotation (Line(points={{-118,-78},{-145,-78}}, color={0,0,127}));
    connect(co2.CO2, PIDco2.u_m) annotation (Line(points={{36.6,70},{64,70},{64,-102},
            {-156,-102},{-156,-90}}, color={0,0,127}));
    connect(PIDco2.y, scale2.u) annotation (Line(points={{-145,-78},{-130,-78},{-130,
            -192},{-112,-192}}, color={0,0,127}));
    connect(PIDt.u_s, tstp)
      annotation (Line(points={{-168,-136},{-226,-136}}, color={0,0,127}));
    connect(PIDco2.u_s, co2stp)
      annotation (Line(points={{-168,-78},{-226,-78}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1)),
      experiment(
        Tolerance=1e-09,
        __Dymola_fixedstepsize=1,
        __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end ZoneCO2R1C1PID;
  annotation (uses(
      Modelica(version="3.2.2"),
      Buildings(version="5.0.1"),
      OU44(version="2")));
end MShootBS2019;
