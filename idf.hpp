#include "simucpp.hpp"
using namespace simucpp;
using namespace zhnmat;
using namespace std;

/*
最小二乘法
输入: 单线 y,x1,x2,x3
输出: 总线 θ
*/
class LeastSquare: public PackModule
{
public:
    LeastSquare() {};
    LeastSquare(Simulator *sim, uint order) : _order(order) {
        Initialize(sim);
    };
    virtual PUnitModule Get_InputPort(int n) const {
        if (n==0) return iny;
        if ((n>=1) && (n<=_order)) return inx[n-1];
        return nullptr;
    };
    virtual PMatModule Get_OutputBus(int n) const {
        if (n==0) return msstheta;
        return nullptr;
    };
    void Set_InitialValue(Mat value) {
        msstheta->Set_InitialValue(value);
    }
    void Set_SampleTime(double time) {
        msstheta->Set_SampleTime(time);
        mssP->Set_SampleTime(time);
    }
    void Initialize(Simulator *sim) {
        SUGain(iny, sim);
        inx = new UGain*[_order];
        for (uint i = 0; i < _order; i++)
            inx[i] = new UGain(sim, "inx");
        msstheta = new MStateSpace(sim, BusSize(_order, 1), SIMUCPP_DISCRETE, "msstheta");
        mssP = new MStateSpace(sim, BusSize(_order, _order), SIMUCPP_DISCRETE, "mssP");
        mxX = new Mux(sim, BusSize(_order, 1), "mxX");
        mxY = new Mux(sim, BusSize(1, 1), "mxY");
        misoK = new MFcnMISO(sim, BusSize(_order, 1), "misoK");
        misoP = new MFcnMISO(sim, BusSize(_order, _order), "misoP");
        misoTheta = new MFcnMISO(sim, BusSize(_order, 1), "misoTheta");
        // mxX := x(k)
        for (uint i = 0; i < _order; i++)
            sim->connectU(inx[i], mxX, BusSize(i, 0));
        // misoK := K(k)
        sim->connectM(mxX, misoK);  // mxX := x(k)
        sim->connectM(mssP, misoK);  // mssP := P(k-1)
        // misoP := P(k)
        sim->connectM(misoK, misoP);  // misoK := K(k)
        sim->connectM(mxX, misoP);  // mxX := x(k)
        sim->connectM(mssP, misoP);  // mssP := P(k-1)
        sim->connectM(misoP, mssP);  // misoP := P(k)
        // misoTheta := θ(k)
        sim->connectU(iny, mxY, BusSize(0, 0));  // iny := y, mxY := y(k)
        sim->connectM(misoK, misoTheta);  // misoK := K(k)
        sim->connectM(mxY, misoTheta);  // mxY := y(k)
        sim->connectM(mxX, misoTheta);  // mxX := x(k)
        sim->connectM(msstheta, misoTheta);  // msstheta := θ(k-1)
        sim->connectM(misoTheta, msstheta);  // msstheta := θ(k-1)

        mssP->Set_InitialValue(eye(_order)*1e6);
        misoK->Set_Function([=](Mat *u){
            Mat ans = u[1] * u[0];  // ans = P(k-1)x(k)
            Mat den = u[0].T() * ans;  // den = x(k)^TP(k-1)x(k)
            double c = 1 / (den.at(0, 0) + _lambda);  // c = 1/(λ+x(k)^TP(k-1)x(k))
            return ans * c;  // K(k) = (1/(1+x(k)^TP(k-1)x(k))) * P(k-1)x(k)
        });
        misoP->Set_Function([=](Mat *u){
            Mat ans = eye(_order) - u[0] * u[1].T();  // ans = I - K(k)x(k)^T
            return ans * u[2];  // P(k) = (I - K(k)x(k)^T) * P(k-1)
        });
        misoTheta->Set_Function([](Mat *u){
            Mat eM = u[1] - u[2].T() * u[3];  // y(k) - x(k)^Tθ(k-1)
            double err = eM.at(0, 0);  // e = y(k) - x(k)^Tθ(k-1)
            Mat ans = u[0] * err;  // ans = K(k) * e
            return ans + u[3];  // θ(k) = K(k) * (y(k) - x(k)^Tθ(k-1)) + θ(k-1)
        });
    }
    double _lambda=0.8;  // 遗忘因子
    uint _order;  // 待辨识参数个数
private:
    UGain *iny=nullptr;  // y
    UGain **inx=nullptr;  // xi
    MStateSpace *msstheta=nullptr;  // θ(k)
    MStateSpace *mssP=nullptr;  // P(k)
    Mux *mxX=nullptr;  // x=[x1; x2; x3]
    Mux *mxY=nullptr;  // y(k)
    MFcnMISO *misoK=nullptr;  // K(k)
    MFcnMISO *misoP=nullptr;  // P(k)
    MFcnMISO *misoTheta=nullptr;  // θ(k)
};

/*
特征模型参数辨识
输入: 单线 y,u
输出: 总线 θ
*/
class ParamIdentifier: public PackModule
{
public:
    ParamIdentifier() {};
    ParamIdentifier(Simulator *sim){
        Initialize(sim);
    };
    virtual PUnitModule Get_InputPort(int n) const {
        if (n==0) return iny;
        if (n==1) return inu;
        return nullptr;
    };
    virtual PMatModule Get_OutputBus(int n) const {
        if (n==0) return ls->Get_OutputBus(0);
        return nullptr;
    };
    void Set_SampleTime(double time) {
        ls->Set_SampleTime(time);
        udy1->Set_SampleTime(time);
        udy2->Set_SampleTime(time);
        udu->Set_SampleTime(time);
    }
    void Initialize(Simulator *sim){
        SUGain(iny, sim);
        SUGain(inu, sim);
        SUUnitDelay(udy1, sim);
        SUUnitDelay(udy2, sim);
        SUUnitDelay(udu, sim);
        ls = new LeastSquare(sim, 3);
        sim->connectU(iny, udy1);
        sim->connectU(udy1, udy2);
        sim->connectU(inu, udu);
        sim->connectU(iny, ls, 0);
        sim->connectU(udy1, ls, 1);
        sim->connectU(udy2, ls, 2);
        sim->connectU(udu, ls, 3);
        ls->Set_InitialValue(Mat(3, 1, vecdble{2, -1, 0}));
    }
private:
    UGain *iny=nullptr;  // y
    UGain *inu=nullptr;  // u
    UUnitDelay *udy1=nullptr;  // y(k-1)
    UUnitDelay *udy2=nullptr;  // y(k-2)
    UUnitDelay *udu=nullptr;  // u(k-1)
    LeastSquare *ls=nullptr;  //
};
