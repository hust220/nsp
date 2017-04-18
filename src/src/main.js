// The Vue build version to load with the `import` command
// (runtime-only or standalone) has been set in webpack.base.conf with an alias.
import Vue from 'vue'
import VueResource from 'vue-resource'
import ElementUI from 'element-ui'
import 'element-ui/lib/theme-default/index.css'
import App from './App'
import VueRouter from 'vue-router'
import Home from './components/Home'
import Doc from './components/Doc'
import API from './components/API'
import Download from './components/Download'

Vue.use(ElementUI)
Vue.use(VueResource)
Vue.use(VueRouter)

const router = new VueRouter({
  routes: [
    { path: '/', component: Home },
    { path: '/home', component: Home },
    { path: '/doc', component: Doc },
    { path: '/doc/:theme', component: Doc },
    { path: '/doc/:theme/:item', component: Doc },
    { path: '/api', component: API },
    { path: '/api/:theme', component: API },
    { path: '/api/:theme/:item', component: API },
    { path: '/download', component: Download },
  ],

  // scrollBehavior(to, from, savedPosition) {
  //   console.log('scrollBehavior')
  //   if (to.hash) {
  //     return {
  //       selector: to.hash
  //     }
  //   }
  // }
})

var bus = new Vue()

/* eslint-disable no-new */
new Vue({
  el: '#app',
  template: '<App/>',
  render: h => h(App),
  components: { App },
  router,
})
