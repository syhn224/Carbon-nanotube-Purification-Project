#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <ctime>
#include <limits>

struct RangeValue {
    double min_val;
    double avg_val;
    double max_val;
};

class NMP_CNT_Mixture {
private:
    // NMP sabitleri
    inline static constexpr double RHO_NMP   = 1028.0;   // kg/m³ (≈ g/L)
    inline static constexpr double ETA_NMP   = 1.65e-3;  // Pa·s
    inline static constexpr double K_NMP     = 0.165;    // W/m·K
    inline static constexpr double CP_NMP    = 1700.0;   // J/kg·K
    inline static constexpr double SIGMA_NMP = 1e-6;     // S/m

    // Maksimum viskozite sınırı cP cinsinden
    inline static constexpr double MAX_VISCOSITY_CP = 60000.0;

    // CNT özellik aralıkları
    double rho_cnt_min, rho_cnt_avg, rho_cnt_max;
    double k_cnt_min,   k_cnt_avg,   k_cnt_max;
    double cp_cnt_min,  cp_cnt_avg,  cp_cnt_max;
    double ar_min,      ar_avg,      ar_max;
    double phi_c_min,   phi_c_avg,   phi_c_max;

    double cnt_gram;
    double cozelti_litre;
    double c_g_l;      

    RangeValue phi_cnt;   // hacimsel kesir
    RangeValue w_cnt;     // kütlesel kesir

    // Elektriksel iletkenlik için sabit
    inline static constexpr double SIGMA_0 = 800.0;

    // ondalık basamak sayısını bul 
    static int decimalPlaces(double val) {
        if (std::abs(val - std::round(val)) < std::numeric_limits<double>::epsilon() * 10.0)
            return 0;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(9) << val;
        std::string s = oss.str();
        auto pos = s.find('.');
        if (pos == std::string::npos) return 0;
        size_t last = s.find_last_not_of('0');
        if (last == std::string::npos || last <= pos) return 0;
        return static_cast<int>(last - pos);
    }

    // Viskozite hesaplama Pa·s cinsinden
    double calcSingleViscosity(double phi, double aspect_ratio) const {
        double p = aspect_ratio;
        // viskozite uzun lifler için
        double intrinsic_visc = (p > 1.0) ? ((p * p) / (15.0 * (std::log(2.0 * p) - 1.5)) + 2.5) : 2.5;
        // Maksimum dolgu
        double phi_m = std::min(0.20, 3.5 / p);

        if (phi >= phi_m) {
            double effective_phi = phi_m * 0.98; 
            double rel_visc = std::pow(1.0 - (effective_phi / phi_m), -2.0);
            double visc_Pa_s = ETA_NMP * rel_visc;
            double visc_cP = visc_Pa_s * 1000.0;
            if (visc_cP > MAX_VISCOSITY_CP)
                return MAX_VISCOSITY_CP / 1000.0;
            return visc_Pa_s;
        }
        // Genel durum
        double rel_visc = std::pow(1.0 - (phi / phi_m), -intrinsic_visc * phi_m);
        double visc_Pa_s = ETA_NMP * rel_visc;
        double visc_cP = visc_Pa_s * 1000.0;
        if (visc_cP > MAX_VISCOSITY_CP)
            return MAX_VISCOSITY_CP / 1000.0;
        return visc_Pa_s;
    }

    double calcSingleThermalConductivity(double phi, double k_p) const {
        double n = 6.0;   // şekil
        double k_f = K_NMP;
        double num = k_p + (n - 1.0) * k_f - (n - 1.0) * phi * (k_f - k_p);
        double den = k_p + (n - 1.0) * k_f + phi * (k_f - k_p);
        return k_f * (num / den);
    }

    double calcSingleElectricalConductivity(double phi, double phi_c) const {
        if (phi < phi_c) return SIGMA_NMP;
        double sigma_m = SIGMA_0 * std::pow(phi - phi_c, 1.8);
        return (sigma_m > SIGMA_NMP) ? sigma_m : SIGMA_NMP;
    }

    //yorum
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

    //gram CNT, litre çözelti, CNT tipi
    NMP_CNT_Mixture(double cnt_g, double cozelti_l, CNTType type = MIXED_COMMERCIAL)
        : cnt_gram(cnt_g), cozelti_litre(cozelti_l)
    {
        if (cnt_g < 0 || cozelti_l <= 0) {
            throw std::invalid_argument("Girdiler pozitif ve 0'dan buyuk olmalidir!");
        }

        c_g_l = cnt_gram / cozelti_litre;

        // CNT tipine göre özellik aralıklarını belirle
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

        // Hacimsel kesir aralığı
        phi_cnt.min_val = c_g_l / rho_cnt_max;
        phi_cnt.avg_val = c_g_l / rho_cnt_avg;
        phi_cnt.max_val = c_g_l / rho_cnt_min;

        // rho_mix = phi*rho_cnt + (1-phi)*RHO_NMP, birimler kg/m³ = g/L
        double rho_mix_min = (phi_cnt.min_val * rho_cnt_min) + ((1.0 - phi_cnt.min_val) * RHO_NMP);
        double rho_mix_avg = (phi_cnt.avg_val * rho_cnt_avg) + ((1.0 - phi_cnt.avg_val) * RHO_NMP);
        double rho_mix_max = (phi_cnt.max_val * rho_cnt_max) + ((1.0 - phi_cnt.max_val) * RHO_NMP);

        // c_g_l [g/L] / rho_mix [g/L] → kütlesel kesir
        w_cnt.min_val = c_g_l / rho_mix_max;
        w_cnt.avg_val = c_g_l / rho_mix_avg;
        w_cnt.max_val = c_g_l / rho_mix_min;

        if (w_cnt.min_val > w_cnt.max_val) std::swap(w_cnt.min_val, w_cnt.max_val);
    }

    // Özellik erişim
    RangeValue getDensity() const {
        return {
            (phi_cnt.min_val * rho_cnt_min) + ((1.0 - phi_cnt.min_val) * RHO_NMP),
            (phi_cnt.avg_val * rho_cnt_avg) + ((1.0 - phi_cnt.avg_val) * RHO_NMP),
            (phi_cnt.max_val * rho_cnt_max) + ((1.0 - phi_cnt.max_val) * RHO_NMP)
        };
    }

    RangeValue getViscosity() const {
        double v1 = calcSingleViscosity(phi_cnt.min_val, ar_min);
        double v2 = calcSingleViscosity(phi_cnt.avg_val, ar_avg);
        double v3 = calcSingleViscosity(phi_cnt.max_val, ar_max);
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
        // En kötü durum -> en büyük perkolasyon eşiği
        // En iyi durum -> en küçük perkolasyon eşiği
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
        
        RangeValue visc_cP = { visc.min_val * 1000.0, visc.avg_val * 1000.0, visc.max_val * 1000.0 };
        printRowStandard("Dinamik Vizkozite",        visc_cP, "mPa*s (cP)");

        RangeValue k_tc = getThermalConductivity();
        printRowStandard("Isil Iletkenlik",          k_tc, "W/m*K");

        RangeValue cp = getSpecificHeat();
        printRowStandard("Ozgul Isi Kapasitesi",     cp, "J/kg*K");

        RangeValue sig = getElectricalConductivity();
        printRowScientific("Elektriksel Iletkenlik",  sig, "S/m");

        std::cout << "================================================================================\n";
        std::cout << "                      DEGERLERIN MADDE/AKISKAN KIYASLAMASI (YORUM)             \n";
        std::cout << "================================================================================\n";
        std::cout << "* Yogunluk (" << std::fixed << std::setprecision(1) << rho.avg_val << " kg/m3)    : " << getDensityAnalogy(rho.avg_val) << "\n";
        std::cout << "* Vizkozite (" << std::fixed << std::setprecision(1) << visc_cP.avg_val << " cP)      : " << getViscosityAnalogy(visc_cP.avg_val) << "\n";
        std::cout << "* Isil Iletk. (" << std::fixed << std::setprecision(3) << k_tc.avg_val << " W/m*K) : " << getThermalCondAnalogy(k_tc.avg_val) << "\n";
        std::cout << "* Ozgul Isi Cp (" << std::fixed << std::setprecision(1) << cp.avg_val << " J/kg*K) : " << getCpAnalogy(cp.avg_val) << "\n";
        std::cout << "* Elek. Iletk. (" << std::scientific << std::setprecision(2) << sig.avg_val << " S/m) : " << getElecCondAnalogy(sig.avg_val) << "\n";
        std::cout << "================================================================================\n";

        if (visc_cP.avg_val > 10000.0) {
            std::cout << "[UYARI] Bu derisimde (" << c_g_l << " g/L) karisim sivi formunu kaybetmis,\n"
                      << "        yogun bir MACUN / BULAMAC (Slurry / Gel) haline gelmistir!\n"
                      << "================================================================================\n";
        }
    }

    
    static void exportRangeCSV(double min_gram, double max_gram, double volume,
                               CNTType type, const std::string& filename)
    {
        if (min_gram < 0 || max_gram < 0 || volume <= 0 || min_gram > max_gram) {
            throw std::invalid_argument("Gecersiz aralik veya hacim!");
        }

        int max_dec = std::max(decimalPlaces(min_gram), decimalPlaces(max_gram));
        double step;
        if (max_dec <= 0) step = 1.0;
        else if (max_dec == 1) step = 0.1;
        else if (max_dec == 2) step = 0.01;
        else step = 0.001;

        int num_steps = static_cast<int>(std::round((max_gram - min_gram) / step)) + 1;

        namespace fs = std::filesystem;
        fs::path output_path = fs::current_path() / "raporlar" / filename;
        fs::create_directories(output_path.parent_path());

        std::ofstream out(output_path);
        if (!out) {
            throw std::runtime_error("CSV dosyasi olusturulamadi: " + output_path.string());
        }

        out << "CNT_g,c_g_L,phi_min,phi_avg,phi_max,"
            << "rho_min_kgm3,rho_avg_kgm3,rho_max_kgm3,"
            << "visc_min_cP,visc_avg_cP,visc_max_cP,"
            << "k_min_WmK,k_avg_WmK,k_max_WmK,"
            << "cp_min_JkgK,cp_avg_JkgK,cp_max_JkgK,"
            << "sigma_min_Sm,sigma_avg_Sm,sigma_max_Sm\n";
        out << std::fixed << std::setprecision(6);

        for (int i = 0; i < num_steps; ++i) {
            double gram = min_gram + i * step;
            if (gram < 0) continue;

            NMP_CNT_Mixture mix(gram, volume, type);
            writeCSVRow(out, gram, mix);
        }

        out.close();
        std::cout << "\nTaramaya basarili: " << num_steps << " adim hesaplandi.\n"
                  << "CSV dosyasi: " << output_path.string() << "\n";
    }

   
    static void exportRangeCSV_WtPercent(double min_wt, double max_wt, double volume,
                                         CNTType type, const std::string& filename)
    {
        if (min_wt < 0.0 || max_wt <= 0.0 || volume <= 0.0 || min_wt > max_wt || max_wt >= 100.0) {
            throw std::invalid_argument("Gecersiz wt% araligi! (%0 - %99.999 arasi olmalidir)");
        }

        int max_dec = std::max(decimalPlaces(min_wt), decimalPlaces(max_wt));
        double step_wt;
        if (max_dec <= 0) step_wt = 1.0;
        else if (max_dec == 1) step_wt = 0.1;
        else if (max_dec == 2) step_wt = 0.01;
        else step_wt = 0.001;

        int num_steps = static_cast<int>(std::round((max_wt - min_wt) / step_wt)) + 1;

        namespace fs = std::filesystem;
        fs::path output_path = fs::current_path() / "raporlar" / filename;
        fs::create_directories(output_path.parent_path());

        std::ofstream out(output_path);
        if (!out) {
            throw std::runtime_error("CSV dosyasi olusturulamadi: " + output_path.string());
        }

        out << "wt_pct,CNT_g,c_g_L,phi_min,phi_avg,phi_max,"
            << "rho_min_kgm3,rho_avg_kgm3,rho_max_kgm3,"
            << "visc_min_cP,visc_avg_cP,visc_max_cP,"
            << "k_min_WmK,k_avg_WmK,k_max_WmK,"
            << "cp_min_JkgK,cp_avg_JkgK,cp_max_JkgK,"
            << "sigma_min_Sm,sigma_avg_Sm,sigma_max_Sm\n";
        out << std::fixed << std::setprecision(6);

        for (int i = 0; i < num_steps; ++i) {
            double wt = min_wt + i * step_wt;
            if (wt <= 0.0 || wt >= 100.0) continue;


            double gram = approximateGramFromWtPercent(wt, volume, type);
            NMP_CNT_Mixture mix(gram, volume, type);
            writeCSVRow(out, wt, gram, mix);   
        }

        out.close();
        std::cout << "\nTaramaya basarili: " << num_steps << " adim hesaplandi (wt% cinsinden).\n"
                  << "CSV dosyasi: " << output_path.string() << "\n";
    }

private:

    static void writeCSVRow(std::ofstream& out, double gram, const NMP_CNT_Mixture& mix) {
        auto phi = mix.phi_cnt;
        auto rho = mix.getDensity();
        auto visc = mix.getViscosity();   
        auto k = mix.getThermalConductivity();
        auto cp = mix.getSpecificHeat();
        auto sig = mix.getElectricalConductivity();

        out << gram << ","
            << mix.c_g_l << ","
            << phi.min_val << "," << phi.avg_val << "," << phi.max_val << ","
            << rho.min_val << "," << rho.avg_val << "," << rho.max_val << ","
            << visc.min_val * 1000.0 << "," << visc.avg_val * 1000.0 << "," << visc.max_val * 1000.0 << ","
            << k.min_val << "," << k.avg_val << "," << k.max_val << ","
            << cp.min_val << "," << cp.avg_val << "," << cp.max_val << ","
            << sig.min_val << "," << sig.avg_val << "," << sig.max_val << "\n";
    }


    static void writeCSVRow(std::ofstream& out, double wt_pct, double gram, const NMP_CNT_Mixture& mix) {
        auto phi = mix.phi_cnt;
        auto rho = mix.getDensity();
        auto visc = mix.getViscosity();
        auto k = mix.getThermalConductivity();
        auto cp = mix.getSpecificHeat();
        auto sig = mix.getElectricalConductivity();

        out << wt_pct << ","
            << gram << ","
            << mix.c_g_l << ","
            << phi.min_val << "," << phi.avg_val << "," << phi.max_val << ","
            << rho.min_val << "," << rho.avg_val << "," << rho.max_val << ","
            << visc.min_val * 1000.0 << "," << visc.avg_val * 1000.0 << "," << visc.max_val * 1000.0 << ","
            << k.min_val << "," << k.avg_val << "," << k.max_val << ","
            << cp.min_val << "," << cp.avg_val << "," << cp.max_val << ","
            << sig.min_val << "," << sig.avg_val << "," << sig.max_val << "\n";
    }


    static double approximateGramFromWtPercent(double wt_pct, double volume_L, CNTType type) {

        double rho_guess = RHO_NMP;  
        double gram = (wt_pct / (100.0 - wt_pct)) * volume_L * rho_guess;

        for (int iter = 0; iter < 5; ++iter) {
            NMP_CNT_Mixture temp_mix(gram, volume_L, type);

            RangeValue rho_range = temp_mix.getDensity();
            double rho_mix_avg = rho_range.avg_val;   

            double new_gram = (wt_pct / (100.0 - wt_pct)) * volume_L * rho_mix_avg;
            if (std::abs(new_gram - gram) < 1e-6) break;
            gram = new_gram;
        }
        return gram;
    }
};

std::string kisa(double d) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3) << d;
    std::string s = ss.str();
    s.erase(s.find_last_not_of('0') + 1, std::string::npos);
    if (s.back() == '.') s.pop_back();
    return s;
}

template<typename T>
bool safeInput(T& value, const std::string& prompt) {
    std::cout << prompt;
    if (!(std::cin >> value)) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        return false;
    }
   
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return true;
}

int main() {
    try {
        std::cout << "==============================================\n";
        std::cout << "   NMP + CNT SUSPANSIYON HESAPLAMA ARACI\n";
        std::cout << "==============================================\n";
        std::cout << "Secenekler:\n";
        std::cout << "  [1] Tek bir karisim hesapla ve raporla\n";
        std::cout << "  [2] Belirli bir CNT araligini tara ve CSV dosyasina kaydet\n";
        std::cout << "Seciminiz (1-2): ";
        int secim;
        if (!safeInput(secim, "")) throw std::invalid_argument("Gecersiz secim!");

        if (secim == 1) {
            double cnt_mg, cozelti_ml;
            int tur_secimi;

            if (!safeInput(cnt_mg, "1. Eklenen CNT Miktarini girin (mg): "))
                throw std::invalid_argument("Gecersiz mg!");
            if (!safeInput(cozelti_ml, "2. Toplam Cozelti Hacmini girin (mL): "))
                throw std::invalid_argument("Gecersiz mL!");

            double cnt_g = cnt_mg / 1000.0;
            double cozelti_l = cozelti_ml / 1000.0;

            std::cout << "\n3. CNT Turunu Secin:\n";
            std::cout << "   [1] Karisik / Ticari CNT (Tavsiye edilen)\n";
            std::cout << "   [2] Saf Cok Duvarli (MWCNT)\n";
            std::cout << "   [3] Saf Tek Duvarli (SWCNT)\n";
            if (!safeInput(tur_secimi, "Seciminiz (1-3) [Varsayilan 1]: ")) {
                tur_secimi = 1;
            }
            if (tur_secimi < 1 || tur_secimi > 3) tur_secimi = 1;

            NMP_CNT_Mixture liquid(cnt_g, cozelti_l, static_cast<NMP_CNT_Mixture::CNTType>(tur_secimi));
            liquid.printReport();
        }
        else if (secim == 2) {
            std::cout << "\nTarama Tipini Secin:\n";
            std::cout << "  [1] Gram cinsinden CNT miktari araligi\n";
            std::cout << "  [2] Agirlikca yuzde (wt%) araligi\n";
            int alt_secim;
            if (!safeInput(alt_secim, "Seciminiz (1-2): ")) throw std::invalid_argument("Gecersiz alt secim!");

            int tur_secimi;
            std::cout << "\nCNT Turunu Secin:\n";
            std::cout << "   [1] Karisik / Ticari CNT\n";
            std::cout << "   [2] Saf Cok Duvarli (MWCNT)\n";
            std::cout << "   [3] Saf Tek Duvarli (SWCNT)\n";
            if (!safeInput(tur_secimi, "Seciminiz (1-3) [Varsayilan 1]: ")) {
                tur_secimi = 1;
            }
            if (tur_secimi < 1 || tur_secimi > 3) tur_secimi = 1;
            NMP_CNT_Mixture::CNTType cnt_type = static_cast<NMP_CNT_Mixture::CNTType>(tur_secimi);

            double vol_ml;
            if (!safeInput(vol_ml, "Cozelti hacmi (mL): ")) throw std::invalid_argument("Gecersiz mL!");
            double vol_l = vol_ml / 1000.0;

            // Zaman damgası
            std::time_t now = std::time(nullptr);
            char time_buf[20];
            std::strftime(time_buf, sizeof(time_buf), "%Y%m%d_%H%M%S", std::localtime(&now));
            std::string timestamp(time_buf);

            std::string dosya_adi;

            if (alt_secim == 1) {
                double min_g, max_g;
                if (!safeInput(min_g, "Minimum CNT miktari (g): ")) throw std::invalid_argument("Gecersiz deger!");
                if (!safeInput(max_g, "Maksimum CNT miktari (g): ")) throw std::invalid_argument("Gecersiz deger!");

                dosya_adi = "gram_CNT" + std::to_string(tur_secimi) + "_"
                          + kisa(min_g) + "-" + kisa(max_g) + "g_"
                          + std::to_string(static_cast<int>(vol_ml)) + "mL_"
                          + timestamp + ".csv";

                NMP_CNT_Mixture::exportRangeCSV(min_g, max_g, vol_l, cnt_type, dosya_adi);
            }
            else if (alt_secim == 2) {
                double min_wt, max_wt;
                if (!safeInput(min_wt, "Minimum agirlikca yuzde (wt%) [orn: 1.5]: ")) throw std::invalid_argument("Gecersiz deger!");
                if (!safeInput(max_wt, "Maksimum agirlikca yuzde (wt%) [orn: 6.0]: ")) throw std::invalid_argument("Gecersiz deger!");

                dosya_adi = "wt_CNT" + std::to_string(tur_secimi) + "_"
                          + kisa(min_wt) + "-" + kisa(max_wt) + "wt_"
                          + std::to_string(static_cast<int>(vol_ml)) + "mL_"
                          + timestamp + ".csv";

                NMP_CNT_Mixture::exportRangeCSV_WtPercent(min_wt, max_wt, vol_l, cnt_type, dosya_adi);
            }
            else {
                std::cout << "Gecersiz alt secim.\n";
            }
        }
        else {
            std::cout << "Gecersiz secim, program sonlandi.\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "\nHata: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
