# МАТЕРІАЛИ ТА МЕТОДИ

> Версія для аудиту джерел. Позначки [PROJECT SOURCE: …] вказують на конкретні файли репозиторію, з яких узято методичне твердження. Позначки [SOURCE NEEDED: …] означають, що відповідне твердження не слід залишати у фінальному рукописі без додаткового підтвердження.

## Дизайн дослідження

Проведено ретроспективне обчислювальне дослідження, у якому за молекулярною структурою прогнозували агреговане значення інгібувальної активності сполук щодо TRPA1 людини. Основною задачею була регресія для сполук із хімічними каркасами, відсутніми в навчальній частині даних. Окремо оцінювали зв’язок похибки прогнозу з розбіжністю опублікованих значень, хімічною віддаленістю тестової сполуки від навчальної вибірки та способом поділу даних.

[PROJECT SOURCE: STATUS.md; docs/paper_plan.md; results/tables/grid_final_metadata_20260801-152155.json; scripts/analyze_h1_h5.py]

## Джерело даних і критерії відбору

Дані отримано з бази ChEMBL версії 37, дата релізу якої у метаданих проєкту зафіксована як 1 травня 2026 року. Відбирали записи для мішені CHEMBL6007 (TRPA1 людини) з типом показника IC50, точною рівністю `standard_relation = "="` та наявним значенням pChEMBL. Заморожена таблиця рівня окремих вимірювань містила 2196 записів; усі збережені стандартизовані значення IC50 були подані в нмоль/л. Документи, з яких походили записи, охоплювали 2010–2025 роки.

[PROJECT SOURCE: data/raw/trpa1_current_api_metadata.json; data/raw/trpa1_current_api_raw.csv] [SOURCE NEEDED: первинна публікація або офіційний опис ChEMBL і визначення поля pChEMBL]

Перехід від 2196 записів до 1645 сполук був агрегацією, а не вилученням біологічних тестів. Після стандартизації структури записи з однаковим InChIKey об’єднували в один molecule-level рядок, а цільову змінну визначали як медіану збережених pChEMBL. Отримана таблиця містила 1645 унікальних сполук, що походили з 97 assays і 55 документів ChEMBL. Для кожної сполуки збережено стандартизований SMILES, InChIKey, медіану, мінімум, максимум і стандартне відхилення pChEMBL, число вимірювань та хімічний каркас Bemis–Murcko. У наборі було 544 унікальні хімічні каркаси.

[PROJECT SOURCE: data/raw/trpa1_current_api_raw.csv; data/processed/trpa1_primary_dataset.csv; data/assays/assay_table.csv; sources/source_registry.csv; results/tables/grid_final_metadata_20260801-152155.json] [SOURCE NEEDED: оригінальне джерело для визначення каркасів Bemis–Murcko]

### Пізніший аудит походження assay-записів

Після завершення основного benchmark усі 97 assays окремо переглянули й класифікували як `primary` (76 assays), `sensitivity` (9 assays) або `excluded` (12 assays). При зіставленні цієї класифікації з 2196 записами початкового benchmark 1953 записи належали до `primary`, 216 — до `sensitivity`, а 27 записів походили з 12 assays, пізніше позначених як `excluded`. Ці 27 записів стосувалися 22 сполук; для 8 сполук вони були єдиним джерелом активності, тоді як для 14 сполук у наборі також існували записи з `primary` або `sensitivity` assays.

[PROJECT SOURCE: data/raw/trpa1_current_api_raw.csv; data/assays/assay_table.csv; results/tables/FINAL_benchmark_assay_role_summary.csv; results/tables/FINAL_assay_role_record_audit.csv]

Оскільки цей аудит було виконано після формування molecule-level набору й завершення основної сітки моделей, результати benchmark та аналізів H1–H5 відображають повний набір 2196 записів / 1645 сполук і не були ретроспективно перераховані після виключення зазначених 27 записів. Аудит використовували для кількісної оцінки обмежень набору, але не як критерій включення до вже виконаного benchmark.

[PROJECT SOURCE: data/processed/trpa1_primary_dataset.csv; results/tables/grid_final_metadata_20260801-152155.json; data/assays/assay_table.csv; STATUS.md]

## Стандартизація структур і формування цільової змінної

Одиницею аналізу була стандартизована молекулярна структура, ідентифікована за InChIKey. Для кожної структури цільову змінну визначали як медіану всіх збережених значень pChEMBL для IC50. Оскільки всі включені записи мали тип показника IC50, у рукописі це агреговане значення далі позначається як pIC50. Такий показник характеризує медіану літературних вимірювань однієї сполуки, а не її активність у єдиному стандартизованому протоколі.

[PROJECT SOURCE: data/raw/trpa1_current_api_raw.csv; data/processed/trpa1_primary_dataset.csv] [SOURCE NEEDED: офіційне джерело ChEMBL, яке підтверджує трактування pChEMBL для точних IC50 як pIC50]

Наявна проєктна документація описує парсинг структур засобами RDKit, видалення солей або вибір основного фрагмента, побудову канонічного ізомеричного SMILES та InChIKey, після чого записи об’єднували за стандартизованою структурою. Водночас у поточній версії репозиторію немає одного канонічного скрипту, який повністю відтворює перехід від таблиці окремих вимірювань до trpa1_primary_dataset.csv.

[PROJECT SOURCE: docs/methods_fact_sheet.md; data/processed/trpa1_primary_dataset.csv] [SOURCE NEEDED: знайти або відновити точний build-script і підтвердити послідовність видалення солей, вибору фрагмента та генерації InChIKey перед поданням статті] [SOURCE NEEDED: офіційна цитата RDKit]

## Молекулярні представлення

Для кожної сполуки використовували чотири способи числового опису молекулярної структури. Першим був бітовий відбиток Morgan ECFP4, сформований у RDKit із радіусом 2 та довжиною 2048 біт без урахування хіральності.

[PROJECT SOURCE: scripts/Grid_Benchmark.py; results/tables/grid_final_metadata_20260801-152155.json] [SOURCE NEEDED: первинна методична публікація для Morgan/ECFP і офіційна цитата RDKit]

Другим представленням був набір із 15 фізико-хімічних дескрипторів RDKit: молекулярна маса, logP, молярна рефракція, площа полярної поверхні, кількість акцепторів і донорів водневого зв’язку, кількість обертових зв’язків, ароматичних, аліфатичних і насичених кілець, загальна кількість кілець, частка sp3-гібридизованих атомів вуглецю, кількість важких атомів, гетероатомів і LabuteASA.

[PROJECT SOURCE: scripts/Grid_Benchmark.py; scripts/chemberta_benchmark.py] [SOURCE NEEDED: офіційна цитата RDKit і, за потреби, першоджерела окремих дескрипторів]

Третім представленням були 384-вимірні вектори ChemBERTa-CLS, отримані замороженою попередньо навченою моделлю DeepChem/ChemBERTa-77M-MTR із першого токена послідовності; для токенізації SMILES застосовували обрізання до максимальної довжини 128 токенів. Параметри цієї моделі не донавчали на наборі TRPA1 у межах основної сітки порівняння.

[PROJECT SOURCE: scripts/chemberta_benchmark.py; scripts/Grid_Benchmark.py; results/tables/grid_final_metadata_20260801-152155.json] [SOURCE NEEDED: перевірена первинна публікація для використаної версії ChemBERTa]

Четвертим представленням були 768-вимірні вектори MolFormer-Mean, сформовані усередненням прихованих представлень попередньо навченої моделі. У замороженому benchmark-файлі збережено контрольну суму NPZ із цими векторами, однак точний ідентифікатор і ревізію вихідної моделі MolFormer ще потрібно підтвердити за середовищем, у якому виконували генерацію векторів.

[PROJECT SOURCE: scripts/Grid_Benchmark.py; results/tables/grid_final_metadata_20260801-152155.json; docs/methods_fact_sheet.md] [SOURCE NEEDED: точна назва й ревізія моделі MolFormer, використаної для embeddings_all.npz] [SOURCE NEEDED: перевірена первинна публікація MolFormer]

Для додаткових аналізів розбіжності значень, хімічної віддаленості та способу поділу даних використовували Morgan fingerprints у поєднанні з Random Forest і XGBoost. Morgan обрали як спільне представлення для обох алгоритмів, оскільки його scaffold-CV якість була близькою до найкращих комбінацій повної сітки, а аналіз хімічної віддаленості безпосередньо ґрунтувався на Tanimoto similarity між Morgan fingerprints.

[PROJECT SOURCE: results/tables/grid_final_results_20260801-152155.csv; results/tables/FINAL_H1_scaffold_performance.csv; scripts/analyze_h1_h5.py]

## Регресійні алгоритми

Кожне з чотирьох молекулярних представлень подавали двом регресійним алгоритмам на однакових розбиттях даних: випадковому лісу (Random Forest) та XGBoost. Для RandomForestRegressor використовували 500 дерев, random_state = 42 і n_jobs = −1; решту параметрів залишали за замовчуванням бібліотеки scikit-learn, а їх точні значення збережено в метаданих запуску.

[PROJECT SOURCE: scripts/Grid_Benchmark.py; results/tables/grid_final_metadata_20260801-152155.json] [SOURCE NEEDED: первинна публікація Random Forest; офіційна цитата scikit-learn]

Для XGBRegressor використовували 500 дерев, максимальну глибину 6, швидкість навчання 0,05, цільову функцію reg:squarederror, метрику rmse, алгоритм побудови дерев hist, random_state = 42 і n_jobs = −1.

[PROJECT SOURCE: scripts/Grid_Benchmark.py; results/tables/grid_final_metadata_20260801-152155.json] [SOURCE NEEDED: первинна публікація XGBoost]

Для контролю результатів також застосовували базову модель, яка в кожному навчальному фолді передбачала середнє значення цільової змінної.

[PROJECT SOURCE: scripts/Grid_Benchmark.py; results/tables/grid_final_results_20260801-152155.csv]

## Поділ даних і основна оцінка моделей

Основну оцінку проводили методом п’ятикратної перехресної перевірки з групуванням за хімічними каркасами Bemis–Murcko. Сполуки з однаковим каркасом завжди належали до одного фолду, тому каркас тестової сполуки не був представлений у навчальній частині відповідного фолду.

[PROJECT SOURCE: scripts/Grid_Benchmark.py; results/tables/grid_final_fold_assignments_20260801-152155.csv; results/tables/grid_final_metadata_20260801-152155.json] [SOURCE NEEDED: первинна публікація Bemis–Murcko; методичне джерело щодо використання scaffold split для оцінки хімічної екстраполяції]

Використовували три варіанти розподілу каркасів між п’ятьма фолдами. Перший варіант формували детермінованим GroupKFold. Для другого і третього варіантів каркаси незалежно перерозподіляли з початковими значеннями генератора 1001 і 1002, балансуючи кількість сполук у фолдах. Ті самі призначення фолдів застосовували до всіх восьми комбінацій «молекулярне представлення × алгоритм».

[PROJECT SOURCE: scripts/Grid_Benchmark.py; results/tables/grid_final_fold_assignments_20260801-152155.csv; results/tables/grid_final_metadata_20260801-152155.json]

Для кожної сполуки в кожному з трьох варіантів розподілу зберігали прогноз, отриманий моделлю, яка не навчалася на цій сполуці та її каркасі (out-of-fold, OOF). У кожному варіанті розподілу метрики обчислювали на сукупності всіх OOF-прогнозів, після чого подавали середнє та стандартне відхилення між трьома варіантами.

[PROJECT SOURCE: scripts/Grid_Benchmark.py; results/tables/grid_final_oof_20260801-152155.csv; results/tables/grid_final_results_20260801-152155.csv]

Основними показниками якості були середньоквадратична похибка (RMSE), коефіцієнт детермінації R² та ранговий коефіцієнт кореляції Спірмена. Середню абсолютну похибку (MAE) додатково використовували у порівнянні випадкового поділу з поділом за каркасами.

[PROJECT SOURCE: scripts/Grid_Benchmark.py; scripts/analyze_h1_h5.py; results/tables/grid_final_results_20260801-152155.csv; results/tables/FINAL_H5_random_vs_scaffold.csv] [SOURCE NEEDED: за вимогами журналу визначити, чи потрібні окремі стандартні статистичні посилання для RMSE, R², MAE та кореляції Спірмена]

## Аналіз розбіжності опублікованих значень

Щоб технічні дублікати одного тесту не вважали незалежними вимірюваннями, записи спочатку об’єднували медіаною всередині кожної пари «сполука × assay». Розбіжність між тестами оцінювали лише для 393 сполук, представлених щонайменше у двох різних assays.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; results/tables/FINAL_H2_variability_vs_error.csv; results/tables/FINAL_H1_H5_METADATA.json]

У суворішому міждокументному аналізі записи спочатку об’єднували медіаною всередині кожної пари «сполука × документ», після чого розбіжність оцінювали для 52 сполук, наявних щонайменше у двох різних документах.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; results/tables/FINAL_H2_variability_vs_error.csv; results/tables/FINAL_H1_H5_METADATA.json]

Основним показником розбіжності був розмах pIC50 між агрегованими значеннями. Як додаткові показники обчислювали вибіркове стандартне відхилення, медіанне абсолютне відхилення та медіану всіх попарних абсолютних різниць. Для кожної моделі абсолютну OOF-похибку сполуки усереднювали між трьома варіантами scaffold-розподілу, після чого перевіряли її монотонний зв’язок із показниками розбіжності за коефіцієнтом Спірмена.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; results/tables/FINAL_H2_variability_vs_error.csv]

Для кожного з чотирьох показників розбіжності 95% довірчі інтервали коефіцієнта Спірмена оцінювали кластерним bootstrap із повторним вибором цілих хімічних каркасів, 1500 ітерацій. Додатково обчислювали часткову рангову кореляцію з контролем медіанного pIC50 та кількості assays або документів.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; results/tables/FINAL_H2_variability_vs_error.csv; results/tables/FINAL_H1_H5_METADATA.json]

Оскільки для кожної моделі в кожній підвибірці одночасно перевіряли чотири показники розбіжності, їхні p-значення коригували методом Холма окремо в межах кожної комбінації «модель × рівень агрегації» (`different_assays` або `different_documents`). Корекція змінювала оцінку статистичної значущості, але не значення коефіцієнта ρ.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; results/tables/FINAL_H2_variability_vs_error.csv] [SOURCE NEEDED: первинне або методичне джерело для поправки Холма]

## Допоміжне порівняння з класифікаційною постановкою

Для оцінки впливу постановки задачі окремо відтворювали класифікаційний підхід «відомий антагоніст TRPA1 проти випадково відібраної decoy-сполуки». У цьому аналізі використовували 2048-бітні Morgan ECFP4, випадковий стратифікований поділ 80/20 з random_state = 34 та десятикратну стратифіковану перехресну перевірку для ROC AUC.

[PROJECT SOURCE: scripts/MihaiExperimentReplication.py; results/tables/FINAL_H3_mihai_replication_summary.csv] [SOURCE NEEDED: повний перевірений бібліографічний запис Mihai et al., DOI 10.3390/ai1020018]

Порівнювали випадковий ліс, метод опорних векторів із радіально-базисним ядром та повнозв’язну нейронну мережу. Цей аналіз відтворював загальну постановку попередньої роботи, але використовував Morgan fingerprints замість її початкового способу молекулярного опису, тому його розглядали як допоміжне порівняння складності задач, а не як дослівне відтворення всіх методичних деталей.

[PROJECT SOURCE: scripts/MihaiExperimentReplication.py] [SOURCE NEEDED: повний текст Mihai et al. для точного зіставлення початкових molecular descriptors і параметрів моделей]

## Аналіз хімічної області застосовності

Для кожної тестової сполуки в кожному scaffold-фолді визначали максимальну коефіцієнтну подібність Танімото між її Morgan fingerprint і fingerprints усіх сполук навчальної частини. Для кожної сполуки максимальну подібність та абсолютну OOF-похибку усереднювали між трьома варіантами scaffold-розподілу.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; results/tables/FINAL_H4_similarity_vs_error.csv] [SOURCE NEEDED: первинне або методичне джерело для коефіцієнта Танімото у хемоінформатиці]

Зв’язок між максимальною подібністю до навчальної вибірки та абсолютною похибкою оцінювали за коефіцієнтом Спірмена. 95% довірчі інтервали отримували scaffold-кластерним bootstrap із 1500 ітераціями; додатково обчислювали часткову рангову кореляцію з контролем медіанного pIC50.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; results/tables/FINAL_H4_similarity_vs_error.csv; results/tables/FINAL_H1_H5_METADATA.json]

## Порівняння випадкового поділу з поділом за каркасами

Для Random Forest і XGBoost на Morgan fingerprints додатково виконували випадкову п’ятикратну перехресну перевірку на тому самому наборі з 1645 сполук. Використовували три незалежні перемішування з `random_state` 1000, 1001 і 1002, зберігаючи ті самі параметри регресійних алгоритмів, що й в основному scaffold-аналізі.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; results/tables/FINAL_H5_random_oof_morgan.csv; results/tables/FINAL_H5_random_vs_scaffold.csv; data/processed/trpa1_primary_dataset.csv]

Для кожного способу поділу розраховували RMSE, MAE, R² і коефіцієнт Спірмена на сукупності OOF-прогнозів. Метрики обчислювали окремо для кожного з трьох варіантів розподілу, після чого усереднювали між ними.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; results/tables/FINAL_H5_random_vs_scaffold.csv]

Статистичну невизначеність різниці між random і scaffold split оцінювали окремо для RMSE, MAE та R². Одиницею повторного вибору був цілий хімічний каркас, а не окрема молекула. Для 95% довірчих інтервалів різниці виконували 5000 bootstrap-ітерацій із повторним вибором 544 каркасів. Двобічне p-значення отримували в 10 000 перестановках, у яких для кожного каркаса випадково міняли місцями помилки, отримані за random і scaffold validation. У кожній ітерації метрики спочатку обчислювали для трьох розподілів окремо, а потім усереднювали, як у головній таблиці H5.

[PROJECT SOURCE: scripts/test_h5_validation_difference.py; results/tables/FINAL_H5_significance.csv]

Для кожної моделі три p-значення, отримані для RMSE, MAE та R², коригували методом Холма. Коефіцієнт Спірмена подано описово й не включено до цього множинного тестування.

[PROJECT SOURCE: scripts/test_h5_validation_difference.py; results/tables/FINAL_H5_significance.csv] [SOURCE NEEDED: методичне джерело для cluster bootstrap і permutation test на рівні хімічних каркасів, якщо редакція вимагатиме зовнішнього статистичного посилання]

## Програмне забезпечення та відтворюваність

Основний benchmark виконано в Python 3.12.13 із NumPy 2.0.2, pandas 2.2.2, scikit-learn 1.6.1, SciPy 1.16.3, XGBoost 3.3.0 і RDKit 2026.03.4. Точні параметри моделей, призначення фолдів, OOF-прогнози, контрольні суми вхідних файлів і таблиці результатів збережено в репозиторії проєкту.

[PROJECT SOURCE: results/tables/grid_final_metadata_20260801-152155.json; results/tables/grid_final_fold_assignments_20260801-152155.csv; results/tables/grid_final_oof_20260801-152155.csv; results/tables/grid_final_results_20260801-152155.csv] [SOURCE NEEDED: офіційні програмні цитати для Python-пакетів, які журнал вимагатиме включити до списку літератури]

Повторний запуск аналізів H1–H5 виконано в Python 3.13.13 із NumPy 2.4.6, pandas 3.0.3, SciPy 1.17.1, scikit-learn 1.8.0, XGBoost 3.2.0 і RDKit 2026.03.2. Скрипт `scripts/analyze_h1_h5.py` сформував таблиці H1–H5, random-split OOF-прогнози та метадані з контрольними сумами. Окремий скрипт `scripts/test_h5_validation_difference.py` сформував cluster-bootstrap і permutation результати для порівняння random та scaffold validation.

[PROJECT SOURCE: scripts/analyze_h1_h5.py; scripts/test_h5_validation_difference.py; results/tables/FINAL_H1_H5_METADATA.json; results/tables/FINAL_H5_significance.csv]

Канонічними вихідними файлами додаткових аналізів є `FINAL_H1_scaffold_performance.csv`, `FINAL_H2_variability_vs_error.csv`, `FINAL_H3_mihai_replication_summary.csv`, `FINAL_H4_similarity_vs_error.csv`, `FINAL_H5_random_vs_scaffold.csv`, `FINAL_H5_random_oof_morgan.csv` і `FINAL_H5_significance.csv`.

[PROJECT SOURCE: results/tables/FINAL_H1_scaffold_performance.csv; results/tables/FINAL_H2_variability_vs_error.csv; results/tables/FINAL_H3_mihai_replication_summary.csv; results/tables/FINAL_H4_similarity_vs_error.csv; results/tables/FINAL_H5_random_vs_scaffold.csv; results/tables/FINAL_H5_random_oof_morgan.csv; results/tables/FINAL_H5_significance.csv]

> [SOURCE NEEDED: скрипт `scripts/test_h5_validation_difference.py` не записує окремий файл із версіями середовища; якщо цей запуск виконувався не в тому самому середовищі, що `scripts/analyze_h1_h5.py`, версії треба зафіксувати окремо.]

## Етичні аспекти

Дослідження було вторинним аналізом відкритих хіміко-біологічних даних і не включало роботу з людьми, персональними даними або лабораторними тваринами.

[PROJECT SOURCE: data/raw/trpa1_current_api_raw.csv; data/processed/trpa1_primary_dataset.csv] [SOURCE NEEDED: підтвердити за правилами журналу формулювання щодо необхідності або відсутності необхідності окремого етичного схвалення для такого обчислювального дослідження]
