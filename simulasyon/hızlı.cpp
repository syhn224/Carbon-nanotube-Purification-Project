#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>
#include <stdexcept>

struct RangeValue {
    double min_val;
    double avg_val;
    double max_val;
};

class NMP_CNT_Mixture {
private:
    // NMP Cozuculer Sabitleri (25 C)
    inline static constexpr double RHO_NMP   = 1028.0;   // kg/m3
    inline static constexpr double ETA_NMP   = 1.65e-3;  // Pa*s (1.65 cP)
    inline static constexpr double K_NMP     = 0.165;    // W/m*K
    inline static constexpr double CP_NMP    = 1700.0;   // J/kg*K
    inline static constexpr double SIGMA_NMP = 1e-6;     // S/m

    // CNT Parametreleri
    double rho_cnt_min, rho_cnt_avg, rho_cnt_max;
    double k_cnt_min,   k_cnt_avg,   k_cnt_max;
    double cp_cnt_min,  cp_cnt_avg,  cp_cnt_max;
    double ar_min,      ar_avg,      ar_max;
    double phi_c_min,   phi_c_avg,   phi_c_max;

    double cnt_gram;
    double cozelti_litre;
    double c_g_l;

    RangeValue phi_cnt;
    RangeValue w_cnt;

    // Gercekci Vizkozite Hesabi
    double calcSingleViscosity(double phi, double aspect_ratio) const {
        double p = aspect_ratio;
        double intrinsic_visc = (p > 1.0) ? ((p * p) / (15.0 * (std::log(2.0 * p) - 1.5)) + 2.5) : 2.5;
        double phi_m = std::min(0.20, 3.5 / p);

        if (phi >= phi_m) {
            double effective_phi = std::min(phi, phi_m * 0.98);
            double rel_visc = std::pow(1.0 - (effective_phi / phi_m), -2.0);
            return ETA_NMP * std::min(rel_visc, 50000.0); // Maksimum 50,000 cP macun siniri
        }

        double rel_visc = std::pow(1.0 - (phi / phi_m), -intrinsic_visc * phi_m);
        return ETA_NMP * std::min(rel_visc, 50000.0);
    }

    double calcSingleThermalConductivity(double phi, double k_p) const {
        double n = 6.0;
        double k_f = K_NMP;
        double num = k_p + (n - 1.0) * k_f - (n - 1.0) * phi * (k_f - k_p);
        double den = k_p + (n - 1.0) * k_f + phi * (k_f - k_p);
        return k_f * (num / den);
    }

    double calcSingleElectricalConductivity(double phi, double phi_c) const {
        if (phi < phi_c) return SIGMA_NMP;
        constexpr double SIGMA_0 = 800.0;
        double sigma_m = SIGMA_0 * std::pow(phi - phi_c, 1.8);
        return (sigma_m > SIGMA_NMP) ? sigma_m : SIGMA_NMP;
    }

    // --- KIYASLAMA VE YORUM FONKSIYONLARI ---

    std::string getCpAnalogy(double val) const {
        if (val < 900) return "Aluminyum / Cam benzeri dusuk isi kapasitesi (~900 J/kg*K)";
        if (val < 1300) return "Zeytinyagi / Makine Yagi benzeri (~1200-1400 J/kg*K)";
        if (val < 1800) return "NMP / Aseton / Alkol grubu organik cozuculer (~1700 J/kg*K)";
        if (val < 2500) return "Etilen Glikol / Antifriz grubu (~2400 J/kg*K)";
        return "Su benzeri yuksek isi kapasiteli sivi (~4184 J/kg*K)";
    }

    std::string getViscosityAnalogy(double val_cp) const {
        if (val_cp < 2.0) return "Suya yakin cok akici sivi (~1 cP)";
        if (val_cp < 10.0) return "Ince Motor Yagi / Sut (~5-10 cP)";
        if (val_cp < 100.0) return "Zeytinyagi / Agir Motor Yagilari (~80 cP)";
        if (val_cp < 1000.0) return "Motor Yagi / Akici Bal (~200-800 cP)";
        if (val_cp < 10000.0) return "Sut Receli / Koyu Bal (~2000-8000 cP)";
        if (val_cp < 30000.0) return "Ketcap / Dis Macunu / Domates Salcasi (~10000-25000 cP)";
        return "Koyu Bulamac / Batarya Katot Macunu / Jel (>30000 cP)";
    }

    std::string getDensityAnalogy(double val) const {
        if (val < 950) return "Sut / Benzin araligi";
        if (val < 1010) return "Saf Su (~1000 kg/m3)";
        if (val < 1050) return "Deniz Suyu / Tam Yagli Sut (~1030 kg/m3)";
        if (val < 1100) return "Kan / Etilen Glikol (~1060 kg/m3)";
        if (val < 1300) return "Grit / Gliserin (~1260 kg/m3)";
        return "Yogun Endustriyel Kimyasal / Slurry";
    }

    std::string getThermalCondAnalogy(double val) const {
        if (val < 0.20) return "NMP / Organik Cozuculer (Isil Yalitimci Sivi)";
        if (val < 0.30) return "Zeytinyagi / Motor Yagi civari (~0.17-0.22 W/m*K)";
        if (val < 0.60) return "Saf Su seviyesi (~0.60 W/m*K)";
        if (val < 2.0) return "Cam / Beton seviyesi";
        return "Yuksek Isil Iletken Nanosivi";
    }

    std::string getElecCondAnalogy(double val) const {
        if (val < 1e-4) return "Yalitkan / Saf Su seviyesi";
        if (val < 1e-1) return "Icme Suyu / Zayif Elektrolit (~10^-2 S/m)";
        if (val < 10.0) return "Deniz Suyu / Yari Iletken Sivi (~4-5 S/m)";
        if (val < 1000.0) return "Yuksek Iletken Karbon Macunu / Iletken Polimer";
        return "Metalik Iletkenlik Duzeyi";
    }

public:
    enum CNTType {
        MIXED_COMMERCIAL = 1,
        MWCNT = 2,
        SWCNT = 3
    };

    NMP_CNT_Mixture(double cnt_g, double cozelti_l, CNTType type = MIXED_COMMERCIAL)
        : cnt_gram(cnt_g), cozelti_litre(cozelti_l)
    {
        if (cnt_g < 0 || cozelti_l <= 0) {
            throw std::invalid_argument("Girdiler pozitif ve 0'dan buyuk olmalidir!");
        }

        c_g_l = cnt_gram / cozelti_litre; // g/L

        switch (type) {
            case MIXED_COMMERCIAL:
                rho_cnt_min = 1400.0; rho_cnt_avg = 1600.0; rho_cnt_max = 1800.0;
                k_cnt_min   = 800.0;  k_cnt_avg   = 1200.0; k_cnt_max   = 1500.0;
                cp_cnt_min  = 600.0;  cp_cnt_avg  = 650.0;  cp_cnt_max  = 700.0;
                ar_min      = 30.0;   ar_avg      = 60.0;   ar_max      = 90.0;
                phi_c_min   = 0.0010; phi_c_avg   = 0.0015; phi_c_max   = 0.0020;
                break;

            case MWCNT:
                rho_cnt_min = 1700.0; rho_cnt_avg = 1800.0; rho_cnt_max = 1900.0;
                k_cnt_min   = 1500.0; k_cnt_avg   = 2000.0; k_cnt_max   = 2500.0;
                cp_cnt_min  = 650.0;  cp_cnt_avg  = 700.0;  cp_cnt_max  = 750.0;
                ar_min      = 80.0;   ar_avg      = 120.0;  ar_max      = 160.0;
                phi_c_min   = 0.0008; phi_c_avg   = 0.0010; phi_c_max   = 0.0012;
                break;

            case SWCNT:
                rho_cnt_min = 1200.0; rho_cnt_avg = 1300.0; rho_cnt_max = 1400.0;
                k_cnt_min   = 2500.0; k_cnt_avg   = 3000.0; k_cnt_max   = 3500.0;
                cp_cnt_min  = 550.0;  cp_cnt_avg  = 600.0;  cp_cnt_max  = 650.0;
                ar_min      = 200.0;  ar_avg      = 300.0;  ar_max      = 400.0;
                phi_c_min   = 0.0003; phi_c_avg   = 0.0005; phi_c_max   = 0.0007;
                break;
        }

        // Hacimsel Kesirler
        phi_cnt.min_val = c_g_l / rho_cnt_max;
        phi_cnt.avg_val = c_g_l / rho_cnt_avg;
        phi_cnt.max_val = c_g_l / rho_cnt_min;

        // Karisim Yogunluklari
        double rho_mix_min = (phi_cnt.min_val * rho_cnt_min) + ((1.0 - phi_cnt.min_val) * RHO_NMP);
        double rho_mix_avg = (phi_cnt.avg_val * rho_cnt_avg) + ((1.0 - phi_cnt.avg_val) * RHO_NMP);
        double rho_mix_max = (phi_cnt.max_val * rho_cnt_max) + ((1.0 - phi_cnt.max_val) * RHO_NMP);

        // Kutlesel Kesirler
        w_cnt.min_val = c_g_l / rho_mix_max;
        w_cnt.avg_val = c_g_l / rho_mix_avg;
        w_cnt.max_val = c_g_l / rho_mix_min;

        if (w_cnt.min_val > w_cnt.max_val) std::swap(w_cnt.min_val, w_cnt.max_val);
    }

    RangeValue getDensity() const {
        return {
            (phi_cnt.min_val * rho_cnt_min) + ((1.0 - phi_cnt.min_val) * RHO_NMP),
            (phi_cnt.avg_val * rho_cnt_avg) + ((1.0 - phi_cnt.avg_val) * RHO_NMP),
            (phi_cnt.max_val * rho_cnt_max) + ((1.0 - phi_cnt.max_val) * RHO_NMP)
        };
    }

    RangeValue getViscosity() const {
        double v1 = calcSingleViscosity(phi_cnt.min_val, ar_min) * 1000.0;
        double v2 = calcSingleViscosity(phi_cnt.avg_val, ar_avg) * 1000.0;
        double v3 = calcSingleViscosity(phi_cnt.max_val, ar_max) * 1000.0;
        return { v1, v2, v3 };
    }

    RangeValue getThermalConductivity() const {
        return {
            calcSingleThermalConductivity(phi_cnt.min_val, k_cnt_min),
            calcSingleThermalConductivity(phi_cnt.avg_val, k_cnt_avg),
            calcSingleThermalConductivity(phi_cnt.max_val, k_cnt_max)
        };
    }

    RangeValue getSpecificHeat() const {
        return {
            (w_cnt.min_val * cp_cnt_min) + ((1.0 - w_cnt.min_val) * CP_NMP),
            (w_cnt.avg_val * cp_cnt_avg) + ((1.0 - w_cnt.avg_val) * CP_NMP),
            (w_cnt.max_val * cp_cnt_max) + ((1.0 - w_cnt.max_val) * CP_NMP)
        };
    }

    RangeValue getElectricalConductivity() const {
        return {
            calcSingleElectricalConductivity(phi_cnt.min_val, phi_c_max),
            calcSingleElectricalConductivity(phi_cnt.avg_val, phi_c_avg),
            calcSingleElectricalConductivity(phi_cnt.max_val, phi_c_min)
        };
    }

    void printReport() const {
        auto printRowStandard = [](const std::string& label, const RangeValue& val, const std::string& unit) {
            std::cout << std::left << std::setw(23) << label << " : "
                      << std::right << std::setw(10) << std::fixed << std::setprecision(4) << val.min_val << "  |  "
                      << std::setw(10) << val.avg_val << "  |  "
                      << std::setw(10) << val.max_val << " " << unit << "\n";
        };

        auto printRowScientific = [](const std::string& label, const RangeValue& val, const std::string& unit) {
            std::cout << std::left << std::setw(23) << label << " : "
                      << std::right << std::setw(10) << std::scientific << std::setprecision(2) << val.min_val << "  |  "
                      << std::setw(10) << val.avg_val << "  |  "
                      << std::setw(10) << val.max_val << " " << unit << "\n";
        };

        std::cout << "\n================================================================================\n";
        std::cout << "                 NMP + CNT SUSPANSIYON OZELLIKLERI (ARALIK ANALIZI)            \n";
        std::cout << "================================================================================\n";
        std::cout << "Girdiler: " << cnt_gram << " gram CNT / " << cozelti_litre << " L Cozelti (Derisim: " << c_g_l << " g/L)\n";
        std::cout << "--------------------------------------------------------------------------------\n";
        std::cout << std::left << std::setw(23) << "OZELLIK / PARAMETRE" << " : "
                  << std::right << std::setw(10) << "MINIMUM" << "  |  "
                  << std::setw(10) << "ORTALAMA" << "  |  "
                  << std::setw(10) << "MAKSIMUM" << " BIRIM\n";
        std::cout << "--------------------------------------------------------------------------------\n";

        printRowStandard("Hacimsel CNT Kesri (vol)", {phi_cnt.min_val * 100, phi_cnt.avg_val * 100, phi_cnt.max_val * 100}, "%");
        printRowStandard("Kutlesel CNT Kesri (wt)",  {w_cnt.min_val * 100, w_cnt.avg_val * 100, w_cnt.max_val * 100}, "%");

        RangeValue rho = getDensity();
        printRowStandard("Karisim Yogunlugu",       rho, "kg/m3");

        RangeValue visc = getViscosity();
        printRowStandard("Dinamik Vizkozite",        visc, "mPa*s (cP)");

        RangeValue k_tc = getThermalConductivity();
        printRowStandard("Isil Iletkenlik",          k_tc, "W/m*K");

        RangeValue cp = getSpecificHeat();
        printRowStandard("Ozgul Isi Kapasitesi",     cp, "J/kg*K");

        RangeValue sig = getElectricalConductivity();
        printRowScientific("Elektriksel Iletkenlik",  sig, "S/m");

        // YORUM VE KIYASLAMA BOLUMU
        std::cout << "================================================================================\n";
        std::cout << "                      DEGERLERIN MADDE/AKISKAN KIYASLAMASI (YORUM)             \n";
        std::cout << "================================================================================\n";
        std::cout << "* Yogunluk (" << std::fixed << std::setprecision(1) << rho.avg_val << " kg/m3)    : " << getDensityAnalogy(rho.avg_val) << "\n";
        std::cout << "* Vizkozite (" << std::fixed << std::setprecision(1) << visc.avg_val << " cP)      : " << getViscosityAnalogy(visc.avg_val) << "\n";
        std::cout << "* Isil Iletk. (" << std::fixed << std::setprecision(3) << k_tc.avg_val << " W/m*K) : " << getThermalCondAnalogy(k_tc.avg_val) << "\n";
        std::cout << "* Ozgul Isi Cp (" << std::fixed << std::setprecision(1) << cp.avg_val << " J/kg*K) : " << getCpAnalogy(cp.avg_val) << "\n";
        std::cout << "* Elek. Iletk. (" << std::scientific << std::setprecision(2) << sig.avg_val << " S/m) : " << getElecCondAnalogy(sig.avg_val) << "\n";
        std::cout << "================================================================================\n";

        if (visc.avg_val > 10000.0) {
            std::cout << "[UYARI] Bu derisimde (" << c_g_l << " g/L) karisim sivi formunu kaybetmis,\n"
                      << "        yogun bir MACUN / BULAMAC (Slurry / Gel) haline gelmistir!\n"
                      << "================================================================================\n";
        }
    }
};

int main() {
    try {
        double cnt_g;
        double cozelti_l;
        int tur_secimi;

        std::cout << "1. Eklenen CNT Miktarini girin (Gram): ";
        if (!(std::cin >> cnt_g)) throw std::invalid_argument("Gecersiz Gram!");

        std::cout << "2. Toplam Cozelti Hacmini girin (Litre): ";
        if (!(std::cin >> cozelti_l)) throw std::invalid_argument("Gecersiz Litre!");

        std::cout << "\n3. CNT Turunu Secin:\n";
        std::cout << "   [1] Karisik / Ticari CNT (Tavsiye edilen)\n";
        std::cout << "   [2] Saf Cok Duvarli (MWCNT)\n";
        std::cout << "   [3] Saf Tek Duvarli (SWCNT)\n";
        std::cout << "Seciminiz (1-3) [Varsayilan 1]: ";

        std::cin >> tur_secimi;
        if (tur_secimi < 1 || tur_secimi > 3) tur_secimi = 1;

        NMP_CNT_Mixture liquid(cnt_g, cozelti_l, static_cast<NMP_CNT_Mixture::CNTType>(tur_secimi));
        liquid.printReport();

    } catch (const std::exception& e) {
        std::cerr << "\nHata: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}