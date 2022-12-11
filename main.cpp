#include <iostream>
#include <fstream>
#include "idf.hpp"

int main()
{
    constexpr double a = -0.5;
    constexpr double b = 1.1;
    constexpr double c = 2.1;
    cout << "start." << endl;
    Simulator sim1(1);
    FUConstant(inu, &sim1);
    FUZOH(zoh1, &sim1);
    FUOutput(out1, &sim1);
    FUOutput(out2, &sim1);
    FUOutput(out3, &sim1);
    auto *idf = new ParamIdentifier(&sim1);
    auto *dmxTheta = new DeMux(&sim1, BusSize(3, 1), "dmxTheta");
    auto tf1 = new TransferFcn(&sim1, vecdble{3, 4}, vecdble{1, 5, 10, 6, 4}, "tf1");
    auto dtf2 = new DiscreteTransferFcn(&sim1, vecdble{0.3}, vecdble{1, -0.97}, "dtf2");

    sim1.connectU(inu, tf1, 0);
    sim1.connectU(tf1, 0, dtf2, 0);
    sim1.connectU(dtf2, 0, zoh1);
    sim1.connectU(zoh1, idf, 0);
    sim1.connectU(inu, idf, 1);
    sim1.connectM(idf, 0, dmxTheta);
    sim1.connectU(dmxTheta, BusSize(0, 0), out1);
    sim1.connectU(dmxTheta, BusSize(1, 0), out2);
    sim1.connectU(dmxTheta, BusSize(2, 0), out3);

    zoh1->Set_SampleTime(0.05);
    idf->Set_SampleTime(0.05);
    dtf2->Set_SampleTime(0.05);
    inu->Set_OutValue(10);
    // out1->Set_EnableStore(false);
    // out2->Set_EnableStore(false);
    // out3->Set_EnableStore(false);
    // sim1.Set_EnableStore(false);
    sim1.Initialize();
    double t = 0;
    // sim1.Simulate_FirstStep();
    // while (t<20) {
    //     t = sim1.Get_t();
    //     sim1.Simulate_OneStep();
    // }
    // cout << "a:  " << out1->Get_OutValue() << "    ";
    // cout << "b:  " << out2->Get_OutValue() << "    ";
    // cout << "c:  " << out3->Get_OutValue() << "    ";
    // cout << endl;
    sim1.Simulate();
    sim1.Plot();
    return 0;
}
