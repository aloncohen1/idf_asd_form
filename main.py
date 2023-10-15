import streamlit as st
import streamlit_analytics

def parse_info(info):
    if info == 'כן':
        return True
    elif info == 'לא':
        return False
    else:
        print(f'cant parse {info}')
        return info

def main():

    home_holder = st.empty()
    with home_holder:

        with st.form("my_form"):
            st.title(' טופס מיון - הפרעת דחק חריפה (ASD)')
            st.subheader("A קריטריון", divider='blue')
            cat_a = st.radio('האם התקיימה חשיפה למוות - ע"י איום / פציעה חמורה / אלימות מינית', ["לא","כן"], index=0)
            st.divider()
            st.subheader("B קריטריון", divider='blue')
            st.info('לאבחנה נדרשת נוכחות של 9 מהסתסמינים הבאים, מכלל חמשת האשוכולות הבאים, ובתנאי שהתסמינים החלו להופיע או להחמיר לאחר האירוע שדווח בקריטריון A')
            st.write('**חודרנות**')
            cat_b_a_1 = st.radio('זכרונות מעיקים וחוזרים לא רצוניים וחודרניים של האירוע', ["לא","כן"])
            cat_b_a_2 = st.radio('סיוטי לילה חוזרים הכוללים תוכן או אפקט הקשור לאירוע', ["לא", "כן"])
            cat_b_a_3 = st.radio('תגובה דיסוציאטיביות בה הנבדק חש כאילו נמצא שוב באירוע', ["לא", "כן"])
            cat_b_a_4 = st.radio('מצוקה פסיכולוגית עזה בחשיפה לגירויים המסמלים את האירוע', ["לא", "כן"])

            cat_b_a_list = [parse_info(i) for i in [cat_b_a_1, cat_b_a_2, cat_b_a_3, cat_b_a_4]]

            st.write('**מצ"ח שלילי**')
            cat_b_b = st.radio('חוסר יכולת מתמשך לחוות רגשות חיוביים (אושר, סיפוק, אהבה)', ["לא","כן"])

            cat_b_b_list = [parse_info(i) for i in [cat_b_b]]

            st.write('**הימנעות**')
            cat_b_c_1 = st.radio('מאמץ להימנע מזיכרונות, מחשבות או רגשות המסמלים את האירוע', ["לא","כן"])
            cat_b_c_2 = st.radio('מאמץ להימנע מגירויים חיצוניים המסמלים את האירוע', ["לא", "כן"])

            cat_b_c_list = [parse_info(i) for i in [cat_b_c_1, cat_b_c_2]]

            st.write('**עוררות**')
            cat_b_d_1 = st.radio('הפרעות שינה', ["לא", "כן"])
            cat_b_d_2 = st.radio('התנהגות רגזינית והתפרצויות זעם', ["לא", "כן"])
            cat_b_d_3 = st.radio('דריכות יתר', ["לא", "כן"])
            cat_b_d_4 = st.radio('קשיים בריכוז', ["לא", "כן"])
            cat_b_d_5 = st.radio('תגובת בהלה מוגזמת', ["לא", "כן"])

            cat_b_d_list = [parse_info(i) for i in [cat_b_d_1, cat_b_d_2, cat_b_d_3, cat_b_d_4, cat_b_d_5]]

            st.write('**דיסוציאציה**')
            cat_b_e_1 = st.radio('תחושה שונה של המציאות - הנבדק עצמו או סביבתו', ["לא", "כן"])
            cat_b_e_2 = st.radio('חוסר יכולת לזכור היבטים חשובים מהאירוע', ["לא", "כן"])

            cat_b_e_list = [parse_info(i) for i in [cat_b_e_1, cat_b_e_2]]
            st.divider()

            st.subheader("C קריטריון", divider='blue')
            st.info("התסמינים מופיעים בדרך כלל בסמוך לאירוע הטראומטי, לאבחנה נדרשת התמדה תסמינית של שלושה ימים לפחות (כלומר 48 ש' ומעלה) ועד חודש ימים")
            cat_c = st.slider("משך התסמינים בימים", 0, 30, 0)
            st.divider()

            st.subheader("D קריטריון", divider='blue')
            cat_d = st.radio('קיימת פגיעה תפקודית במישורים משפחתיים / חברתיים / תעסוקתיים / אחר', ["לא", "כן"])
            st.divider()
            st.subheader("E קריטריון", divider='blue')
            cat_e = st.radio('אין מקור ההפרעה מיוחס להשפעות של חומרים (אלכוהול / סמים / תרופות) או למצב רפואי אחר (למשל פגיעת ראש) וכן אין חשד לייתכנות של מצב פסיכוטי',
                     ["לא מיוחס לחומרים", "מיוחס לחומרים"])
            st.divider()

            submitted = st.form_submit_button("סיום")

        if submitted:

            has_asd = "החייל תקין"

            if parse_info(cat_a):
                if any(cat_b_a_list) and any(cat_b_b_list) and any(cat_b_c_list) and any(cat_b_d_list) and any(cat_b_e_list):
                    if sum(cat_b_a_list+cat_b_b_list+cat_b_c_list+cat_b_d_list+cat_b_e_list) >= 9:
                        if cat_c >= 3:
                            if parse_info(cat_d):
                                if cat_e == "לא מיוחס לחומרים":
                                    has_asd = "דחוף לאברבנל"

            home_holder.empty()
            st.title(has_asd)


if __name__ == "__main__":
        streamlit_analytics.start_tracking()
        main()
        streamlit_analytics.stop_tracking(unsafe_password='asd_form')