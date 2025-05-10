import org.opensim.modeling.*;
model=Model("OrthoModel_2kgWeight_scaled.osim");
GH=model.getMarkerSet().get("GH").getLocationInGround(model.initSystem()).getAsMat();
EL=model.getMarkerSet().get("EL").getLocationInGround(model.initSystem()).getAsMat();
EM=model.getMarkerSet().get("EM").getLocationInGround(model.initSystem()).getAsMat();
US2=model.getMarkerSet().get("US2").getLocationInGround(model.initSystem()).getAsMat();
RS2=model.getMarkerSet().get("RS2").getLocationInGround(model.initSystem()).getAsMat();

elbow_cent=(EL+EM)/2;
wrist_cent=(US2+RS2)/2;

L=norm(elbow_cent-GH)+norm(wrist_cent-elbow_cent);

Error=0.02*L*100;

%%
import org.opensim.modeling.*;
model=Model("OrthoModel_2kgWeight_scaled.osim");
GCC=model.getMarkerSet().get("Glenoid_Center").getLocationInGround(model.initSystem()).getAsMat();
GCE=model.getMarkerSet().get("Glenoid_Edge").getLocationInGround(model.initSystem()).getAsMat();
norm(GCC-GCE)
%%
import org.opensim.modeling.*;
adapt=TRCFileAdapter();
adapt.read("transformed_MSabdkg1.trc")
dataTable=adapt.getDataTable(adapt.read("transformed_MSabdkg1.trc"),"markers")

osimTableToStruct(dataTable)
%%

T=osimTableToStruct(Tabel);
T1=osimTableToStruct(Tabel);