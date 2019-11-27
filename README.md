## Задание из курса "Управление вычислениями" на платформе Stepik

Я решаю здесь простую задачу из биоинформатики, сама задача не имеет значения для проекта, это могла быть также задача из эконометрики или другой области, вычислительно сложная

### Технологии
Snakemake, Docker,


1. Snakemake - язык описания нагруженных вычислительных задач, хорошо документирован, использует Python для связывания кода

2. Имеет "низкий порог входа" - интутивно понятен, поэтому легко освоить

3. Удобно отлаживать - не пересчитывает уже обработанные шаги, которые отработали без ошибок

4. На мой взгляд, удобно использовать в проектах, разрабатываемых с нуля, когда можно договориться и передаваемых форматах данных

5. между своими процессами передает только файлы с данными, которыми сам управляет

6. Однако не позволяет использовать многие утилиты сторонних разработчиков, т.к. часть утилит не возвращает файлы, это могут быть запросы к базам данных, файлы в жестко зашитой директории (которыми нет возможности управлять в процессе исполнения), структуры данных JSON,  строковые, и т.п.

7. Нет встроенных средств контейнеризации workflow, а в нагруженных вычислениях это обязательный шаг практически всегда
